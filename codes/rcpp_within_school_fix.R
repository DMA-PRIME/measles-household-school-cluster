# ==============================================================================
# CORRECTED: Within-School Transmission (density-dependent contact sampling)
# ==============================================================================
#
# FIX: Contacts are now sampled from ALL present students (any disease state),
# not just S/V. This correctly models survey-derived contact rates which
# represent total classroom encounters. Contacts that land on non-susceptible
# individuals (E, P, Ra, R) are "wasted" — no transmission occurs.
#
# CHANGE SUMMARY:
#   - Added class_members / all_members pools (full present population)
#   - Contact sampling draws from class_members / all_members
#   - Transmission only occurs if sampled target is S or V
#   - Contact tracing records ALL contacts (regardless of target state)
#
# Replace the cpp_school_transmission_contacts cppFunction block in
# rcpp_transmission.R with this corrected version.
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

  // -----------------------------------------------------------------------
  // FIX: Build TWO sets of pools:
  //   1. class_members / all_members: ALL present students (for sampling who
  //      the infector contacts — represents survey-derived contact rates)
  //   2. State checking at transmission time: only S/V can become infected
  // -----------------------------------------------------------------------

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
