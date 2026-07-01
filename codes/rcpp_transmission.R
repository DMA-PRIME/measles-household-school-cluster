# ==============================================================================
# RCPP TRANSMISSION FUNCTIONS
# ==============================================================================
# File: rcpp_transmission.R
# Contains: Compiled C++ functions for fast transmission calculations
# Note: Using rebuild=TRUE to avoid cached version conflicts in parallel workers
# ==============================================================================

library(Rcpp)

# Force rebuild to avoid cache issues with parallel workers
# This adds ~2-3 seconds to loading but ensures consistent function signatures
RCPP_REBUILD <- TRUE

# ==============================================================================
# Within-school transmission (contact-based with tracking)
# ==============================================================================



cppFunction('

List cpp_school_transmission_contacts(

    IntegerVector student_id,

    IntegerVector class_id,

    CharacterVector state,

    LogicalVector is_vaccinated,

    LogicalVector is_isolated,

    LogicalVector is_quarantined,

    double c_within,

    double c_between,

    double p_within,

    double p_between,

    double prodromal_mult,

    double rash_mult,

    double vaccine_reduction,

    double vaccine_efficacy

) {

  int n = student_id.size();

  std::vector<int> new_exposures;

  std::vector<int> breakthrough_cases;

  std::vector<int> contact_infector_ids;

  std::vector<int> contact_target_ids;



  // Identify infectious individuals (P, QP, Ra — not isolated/quarantined)

  std::vector<int> prodromal_idx;

  std::vector<int> rash_idx;



  for (int i = 0; i < n; i++) {

    if (!is_isolated[i] && !is_quarantined[i]) {

      std::string st = Rcpp::as<std::string>(state[i]);

      if (st == "P" || st == "QP") {

        prodromal_idx.push_back(i);

      } else if (st == "Ra") {

        rash_idx.push_back(i);

      }

    }

  }



  if (prodromal_idx.empty() && rash_idx.empty()) {

    return List::create(

      Named("new_exposures") = IntegerVector(0),

      Named("breakthrough_cases") = IntegerVector(0),

      Named("contact_infector_ids") = IntegerVector(0),

      Named("contact_target_ids") = IntegerVector(0)

    );

  }



  // Full present population by class (for contact sampling)

  std::map<int, std::vector<int>> class_members;

  std::vector<int> all_members;



  for (int i = 0; i < n; i++) {

    if (!is_isolated[i] && !is_quarantined[i]) {

      class_members[class_id[i]].push_back(i);

      all_members.push_back(i);

    }

  }



  if (all_members.empty()) {

    return List::create(

      Named("new_exposures") = IntegerVector(0),

      Named("breakthrough_cases") = IntegerVector(0),

      Named("contact_infector_ids") = IntegerVector(0),

      Named("contact_target_ids") = IntegerVector(0)

    );

  }



  // Pre-compute which individuals are susceptible (S or V)

  std::vector<bool> is_susceptible(n, false);

  for (int i = 0; i < n; i++) {

    std::string st = Rcpp::as<std::string>(state[i]);

    if ((st == "S" || st == "V") && !is_isolated[i] && !is_quarantined[i]) {

      is_susceptible[i] = true;

    }

  }



  std::set<int> newly_exposed;

  std::set<int> new_breakthrough;



  // Process infectious individuals

  auto process_infector = [&](int inf_idx, double infectiousness_mult, double p_w, double p_b) {

    int inf_class = class_id[inf_idx];

    int inf_id = student_id[inf_idx];



    // Within-class contacts: sample from ALL classmates (not just susceptibles)

    int n_contacts_within = R::rpois(c_within);

    std::vector<int>& same_class = class_members[inf_class];



    if (!same_class.empty() && n_contacts_within > 0) {

      for (int c = 0; c < n_contacts_within; c++) {

        int target_idx = same_class[rand() % same_class.size()];



        // Skip self-contact

        if (target_idx == inf_idx) continue;



        // Record contact for tracing (all contacts, not just susceptible)

        contact_infector_ids.push_back(inf_id);

        contact_target_ids.push_back(student_id[target_idx]);



        // Transmission only possible if target is susceptible (S or V)

        if (!is_susceptible[target_idx]) continue;

        if (newly_exposed.find(target_idx) != newly_exposed.end()) continue;



        double trans_prob = p_w * infectiousness_mult;



        std::string target_state = Rcpp::as<std::string>(state[target_idx]);

        if (target_state == "V") {

          if (R::runif(0, 1) < vaccine_efficacy) continue;

          trans_prob *= (1.0 - vaccine_reduction);

        }



        if (R::runif(0, 1) < trans_prob) {

          newly_exposed.insert(target_idx);

          if (is_vaccinated[target_idx]) {

            new_breakthrough.insert(target_idx);

          }

        }

      }

    }



    // Between-class contacts: sample from ALL students in school

    int n_contacts_between = R::rpois(c_between);

    if (!all_members.empty() && n_contacts_between > 0) {

      for (int c = 0; c < n_contacts_between; c++) {

        int target_idx = all_members[rand() % all_members.size()];



        // Skip same-class and self

        if (class_id[target_idx] == inf_class) continue;

        if (target_idx == inf_idx) continue;



        // Record contact for tracing

        contact_infector_ids.push_back(inf_id);

        contact_target_ids.push_back(student_id[target_idx]);



        // Transmission only if target is susceptible

        if (!is_susceptible[target_idx]) continue;

        if (newly_exposed.find(target_idx) != newly_exposed.end()) continue;



        double trans_prob = p_b * infectiousness_mult;



        std::string target_state = Rcpp::as<std::string>(state[target_idx]);

        if (target_state == "V") {

          if (R::runif(0, 1) < vaccine_efficacy) continue;

          trans_prob *= (1.0 - vaccine_reduction);

        }



        if (R::runif(0, 1) < trans_prob) {

          newly_exposed.insert(target_idx);

          if (is_vaccinated[target_idx]) {

            new_breakthrough.insert(target_idx);

          }

        }

      }

    }

  };



  for (int i : prodromal_idx) {

    process_infector(i, prodromal_mult, p_within, p_between);

  }



  for (int i : rash_idx) {

    process_infector(i, rash_mult, p_within, p_between);

  }



  return List::create(

    Named("new_exposures") = IntegerVector(newly_exposed.begin(), newly_exposed.end()),

    Named("breakthrough_cases") = IntegerVector(new_breakthrough.begin(), new_breakthrough.end()),

    Named("contact_infector_ids") = wrap(contact_infector_ids),

    Named("contact_target_ids") = wrap(contact_target_ids)

  );

}

', rebuild = RCPP_REBUILD)

# ==============================================================================
# Between-school transmission
# ==============================================================================

cppFunction('
List cpp_between_school_transmission(
    IntegerVector infector_ids,
    IntegerVector infector_school_ids,
    CharacterVector infector_states,
    LogicalVector infector_vaccinated,
    List target_pools,
    IntegerVector connected_schools,
    NumericVector edge_weights,
    double c_between_base,
    double p_base,
    double prodromal_mult,
    double rash_mult,
    double vaccine_reduction,
    double vaccine_efficacy
) {
  int n_infectors = infector_ids.size();
  int n_target_schools = target_pools.size();
  
  std::vector<int> exposed_student_ids;
  std::vector<int> exposed_school_ids;
  std::vector<bool> is_breakthrough;
  
  if (n_infectors == 0 || n_target_schools == 0) {
    return List::create(
      Named("exposed_student_ids") = IntegerVector(0),
      Named("exposed_school_ids") = IntegerVector(0),
      Named("is_breakthrough") = LogicalVector(0)
    );
  }
  
  std::set<int> already_exposed;
  
  for (int i = 0; i < n_infectors; i++) {
    std::string inf_state = Rcpp::as<std::string>(infector_states[i]);
    double infectiousness_mult = (inf_state == "P") ? prodromal_mult : rash_mult;
    
    for (int t = 0; t < n_target_schools; t++) {
      double weight = edge_weights[t];
      double effective_contact_rate = c_between_base * weight;
      
      int n_contacts = R::rpois(effective_contact_rate);
      if (n_contacts == 0) continue;
      
      DataFrame target_df = as<DataFrame>(target_pools[t]);
      IntegerVector target_ids = target_df["student_id"];
      IntegerVector target_school_ids = target_df["school_id"];
      CharacterVector target_states = target_df["state"];
      LogicalVector target_vaccinated = target_df["is_vaccinated"];
      
      int n_targets = target_ids.size();
      if (n_targets == 0) continue;
      
      for (int c = 0; c < n_contacts; c++) {
        int target_idx = rand() % n_targets;
        int target_id = target_ids[target_idx];
        
        if (already_exposed.find(target_id) != already_exposed.end()) continue;
        
        double trans_prob = p_base * infectiousness_mult;
        
        std::string target_state = Rcpp::as<std::string>(target_states[target_idx]);
        bool is_vacc = target_vaccinated[target_idx];
        
        if (target_state == "V") {
          if (R::runif(0, 1) < vaccine_efficacy) continue;
          trans_prob *= (1.0 - vaccine_reduction);
        }
        
        if (R::runif(0, 1) < trans_prob) {
          already_exposed.insert(target_id);
          exposed_student_ids.push_back(target_id);
          exposed_school_ids.push_back(target_school_ids[target_idx]);
          is_breakthrough.push_back(is_vacc);
        }
      }
    }
  }
  
  return List::create(
    Named("exposed_student_ids") = wrap(exposed_student_ids),
    Named("exposed_school_ids") = wrap(exposed_school_ids),
    Named("is_breakthrough") = wrap(is_breakthrough)
  );
}
', rebuild = RCPP_REBUILD)


# ==============================================================================
# Quarantine application with contact history
# ==============================================================================

cppFunction('
List cpp_apply_quarantine_with_history(
    IntegerVector student_id,
    CharacterVector state,
    LogicalVector is_quarantined,
    LogicalVector is_vaccinated,
    LogicalVector newly_isolated,
    IntegerVector contact_history_infector,
    IntegerVector contact_history_target,
    double quarantine_efficacy
) {
  int n = student_id.size();
  
  std::set<int> isolated_ids;
  for (int i = 0; i < n; i++) {
    if (newly_isolated[i]) {
      isolated_ids.insert(student_id[i]);
    }
  }
  
  if (isolated_ids.empty()) {
    return List::create(
      Named("quarantine_ids") = IntegerVector(0),
      Named("quarantine_states") = CharacterVector(0)
    );
  }
  
  std::set<int> contacts_of_isolated;
  int n_contacts = contact_history_infector.size();
  
  for (int i = 0; i < n_contacts; i++) {
    if (isolated_ids.find(contact_history_infector[i]) != isolated_ids.end()) {
      contacts_of_isolated.insert(contact_history_target[i]);
    }
  }
  
  std::map<int, int> id_to_idx;
  for (int i = 0; i < n; i++) {
    id_to_idx[student_id[i]] = i;
  }
  
  std::vector<int> quarantine_ids;
  std::vector<std::string> quarantine_states;
  
  for (int contact_id : contacts_of_isolated) {
    if (id_to_idx.find(contact_id) == id_to_idx.end()) continue;
    
    int idx = id_to_idx[contact_id];
    
    if (is_quarantined[idx]) continue;
    
    std::string current_state = Rcpp::as<std::string>(state[idx]);
    
    if (current_state != "S" && current_state != "V" && 
        current_state != "E" && current_state != "P") continue;
    
    if (R::runif(0, 1) > quarantine_efficacy) continue;
    
    std::string new_state;
    if (current_state == "S") {
      new_state = "QS";
    } else if (current_state == "V") {
      new_state = "QS";
    } else if (current_state == "E") {
      new_state = "QE";
    } else if (current_state == "P") {
      new_state = "QP";
    }
    
    quarantine_ids.push_back(contact_id);
    quarantine_states.push_back(new_state);
  }
  
  return List::create(
    Named("quarantine_ids") = wrap(quarantine_ids),
    Named("quarantine_states") = wrap(quarantine_states)
  );
}
', rebuild = RCPP_REBUILD)


# ==============================================================================
# Household transmission (FAST C++ version)
# ==============================================================================

cppFunction('
List cpp_household_transmission(
    NumericVector student_hh_id,
    CharacterVector student_state,
    LogicalVector student_is_vaccinated,
    IntegerVector student_school_idx,
    IntegerVector student_local_idx,
    NumericVector hh_member_hh_id,
    IntegerVector hh_member_id,
    CharacterVector hh_member_state,
    LogicalVector hh_member_is_vaccinated,
    LogicalVector hh_member_is_student,
    double hh_transmission_prob,
    double vaccine_efficacy,
    double vaccine_infectiousness_reduction
) {
  int n_students = student_hh_id.size();
  int n_hh_members = hh_member_hh_id.size();
  
  // Result vectors
  std::vector<int> exposed_student_school_idx;
  std::vector<int> exposed_student_local_idx;
  std::vector<bool> exposed_student_breakthrough;
  
  std::vector<int> exposed_hh_member_id;
  std::vector<bool> exposed_hh_member_breakthrough;
  
  // Build map of household -> infectious count
  std::map<double, int> hh_infectious_count;
  
  // Count infectious students per household
  for (int i = 0; i < n_students; i++) {
    std::string st = Rcpp::as<std::string>(student_state[i]);
    if (st == "P" || st == "Ra" || st == "Iso" || st == "QP") {
      double hh = student_hh_id[i];
      if (!ISNA(hh)) {
        hh_infectious_count[hh]++;
      }
    }
  }
  
  // Count infectious non-student household members
  for (int i = 0; i < n_hh_members; i++) {
    if (hh_member_is_student[i]) continue;  // Skip students (already counted)
    
    std::string st = Rcpp::as<std::string>(hh_member_state[i]);
    if (st == "P" || st == "Ra" || st == "Iso" || st == "QP") {
      double hh = hh_member_hh_id[i];
      if (!ISNA(hh)) {
        hh_infectious_count[hh]++;
      }
    }
  }
  
  if (hh_infectious_count.empty()) {
    return List::create(
      Named("exposed_student_school_idx") = IntegerVector(0),
      Named("exposed_student_local_idx") = IntegerVector(0),
      Named("exposed_student_breakthrough") = LogicalVector(0),
      Named("exposed_hh_member_id") = IntegerVector(0),
      Named("exposed_hh_member_breakthrough") = LogicalVector(0),
      Named("n_infectious_households") = 0
    );
  }
  
  // Process susceptible students
  for (int i = 0; i < n_students; i++) {
    double hh = student_hh_id[i];
    if (ISNA(hh)) continue;
    
    // Check if household has infectious members
    auto it = hh_infectious_count.find(hh);
    if (it == hh_infectious_count.end()) continue;
    
    int n_infectious = it->second;
    
    std::string st = Rcpp::as<std::string>(student_state[i]);
    if (st != "S" && st != "V") continue;
    
    // Determine effective transmission probability
    double p_eff = hh_transmission_prob;
    bool is_vacc = student_is_vaccinated[i];
    
    if (st == "V") {
      // Vaccine protection check
      if (R::runif(0, 1) < vaccine_efficacy) {
        continue;  // Protected by vaccine
      }
      p_eff = p_eff * (1.0 - vaccine_infectiousness_reduction);
    }
    
    // Mass action: P(infection) = 1 - (1-p)^n_infectious
    double p_infection = 1.0 - pow(1.0 - p_eff, n_infectious);
    
    if (R::runif(0, 1) < p_infection) {
      exposed_student_school_idx.push_back(student_school_idx[i]);
      exposed_student_local_idx.push_back(student_local_idx[i]);
      exposed_student_breakthrough.push_back(is_vacc);
    }
  }
  
  // Process susceptible non-student household members
  for (int i = 0; i < n_hh_members; i++) {
    if (hh_member_is_student[i]) continue;
    
    double hh = hh_member_hh_id[i];
    if (ISNA(hh)) continue;
    
    auto it = hh_infectious_count.find(hh);
    if (it == hh_infectious_count.end()) continue;
    
    int n_infectious = it->second;
    
    std::string st = Rcpp::as<std::string>(hh_member_state[i]);
    if (st != "S" && st != "V") continue;
    
    double p_eff = hh_transmission_prob;
    bool is_vacc = hh_member_is_vaccinated[i];
    
    if (st == "V") {
      if (R::runif(0, 1) < vaccine_efficacy) {
        continue;
      }
      p_eff = p_eff * (1.0 - vaccine_infectiousness_reduction);
    }
    
    double p_infection = 1.0 - pow(1.0 - p_eff, n_infectious);
    
    if (R::runif(0, 1) < p_infection) {
      exposed_hh_member_id.push_back(hh_member_id[i]);
      exposed_hh_member_breakthrough.push_back(is_vacc);
    }
  }
  
  return List::create(
    Named("exposed_student_school_idx") = wrap(exposed_student_school_idx),
    Named("exposed_student_local_idx") = wrap(exposed_student_local_idx),
    Named("exposed_student_breakthrough") = wrap(exposed_student_breakthrough),
    Named("exposed_hh_member_id") = wrap(exposed_hh_member_id),
    Named("exposed_hh_member_breakthrough") = wrap(exposed_hh_member_breakthrough),
    Named("n_infectious_households") = (int)hh_infectious_count.size()
  );
}
', rebuild = RCPP_REBUILD)


# ==============================================================================
# Update household member disease states (FAST C++ version)
# ==============================================================================

cppFunction('
List cpp_update_hh_disease_states(
    CharacterVector state,
    LogicalVector is_student,
    IntegerVector time_in_state,
    IntegerVector latent_duration,
    IntegerVector prodromal_duration,
    IntegerVector rash_duration,
    IntegerVector time_since_prodromal
) {
  int n = state.size();
  
  CharacterVector new_state = clone(state);
  IntegerVector new_time_in_state = clone(time_in_state);
  IntegerVector new_time_since_prodromal = clone(time_since_prodromal);
  
  for (int i = 0; i < n; i++) {
    if (is_student[i]) continue;  // Skip students
    
    std::string st = Rcpp::as<std::string>(state[i]);
    
    // E -> P transition
    if (st == "E" && time_in_state[i] >= latent_duration[i]) {
      new_state[i] = "P";
      new_time_in_state[i] = 0;
      new_time_since_prodromal[i] = 0;
    }
    // P -> Ra transition
    else if (st == "P" && time_in_state[i] >= prodromal_duration[i]) {
      new_state[i] = "Ra";
      new_time_in_state[i] = 0;
    }
    // Ra -> R transition
    else if (st == "Ra" && time_in_state[i] >= rash_duration[i]) {
      new_state[i] = "R";
      new_time_in_state[i] = 0;
    }
    else {
      // Increment time in state
      new_time_in_state[i] = time_in_state[i] + 1;
      
      // Track time since prodromal for infectious individuals
      if (st == "P" || st == "Ra") {
        if (IntegerVector::is_na(time_since_prodromal[i])) {
          new_time_since_prodromal[i] = 1;
        } else {
          new_time_since_prodromal[i] = time_since_prodromal[i] + 1;
        }
      }
    }
  }
  
  return List::create(
    Named("state") = new_state,
    Named("time_in_state") = new_time_in_state,
    Named("time_since_prodromal") = new_time_since_prodromal
  );
}
', rebuild = RCPP_REBUILD)


cat("Rcpp transmission functions compiled successfully.\n")
cat("  - cpp_school_transmission_contacts: School-based transmission\n")
cat("  - cpp_household_transmission: Household-based transmission (NEW)\n")
cat("  - cpp_update_hh_disease_states: Update household member states (NEW)\n")
