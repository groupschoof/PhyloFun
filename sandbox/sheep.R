library(gRain);

# Constants
sheep.node_types <- c(~duplication, ~speciation)
sheep.s_speciation <- 2
sheep.s_duplication <- 3

# Exponent of Molecular-Function-Transition-Probability
# dependent on phylogenetic Distance b.
# b is defined as the most parsimonious number of AA-Changes
# along a given branch.
sheep.d <- function(b) {
  1.5 - ( 1.0 / (1.0 + exp(-b)));
}

# Calculates the probability of 
# molecular-function-transition 
# from some molecular function m to n.
# This probability depends only on the 
# phylogenetic type of the tree-node t_n, where
# this probability is to be calculated, and
# on the distance l_g, number of edges, between 
# the two molecular functions in the Gene-Ontology
# Directed Acyclic Graph (GO-DAG).
# For self-transition (m equals n) 
# q is (1 / r) ^ (s / 2), where
# r is the cardinality of the event-space,
# the number of assignable GO-Terms.
sheep.q <- function(l_g, t_n=~speciation, r = 0) {
  s <- if(t_n == ~speciation){ sheep.s_speciation }else{ sheep.s_duplication }

  if(l_g > 0){
    l_g ^ (-s)
  }else{
    (1.0 / r) ^ (s / 2)
  }
}

# child_state_index is needed to know when to infer 
# probability of self-transition and when not.
# child_state_index is the molecular function m we calculate 
# the joint probability of transition to.
sheep.transition_prob <- function(ancestor_states, 
    child_state_index, 
    states_distances, 
    aa_changes_branch, 
    ancestor_node_type = ~speciation) 
{
  # Clone ancestor's random-vector
  transition_probs <- ancestor_states
  for(i in 1:length(ancestor_states)) {
    transition_probs[i] <- sheep.q(states_distances[i, child_state_index],
        ancestor_node_type, length(ancestor_states))
  }
  return(transition_probs)
}


# Calculates the state-probabilities, 
# as a conditional probability table, 
# given parent node's event-space.
# To enable usage with gRain, the table is
# given as a single vector. 
# gRain knows, at what index to break it into rows.
sheep.child_cond_prob_matrix <- function(ancestor_states, 
    states_distances, 
    aa_changes_branch, 
    ancestor_node_type = ~speciation
    ) 
{
  child_states = vector(mode="numeric")
  for(i in 1:length(ancestor_states)) {
    child_states <- append(child_states, 
        sheep.transition_prob(ancestor_states, i, 
            states_distances, aa_changes_branch, 
            ancestor_node_type)
        ) 
  }
  return(child_states)
}
