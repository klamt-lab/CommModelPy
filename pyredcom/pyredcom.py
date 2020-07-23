#!/usr/bin/env python3
#
# Copyright 2020 PSB
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""
This module contains all major class and function definitions of
pyredcom. pyredcom's source code is intersected in two parts:

1. The definition of all dataclasses which describe single models and a community.
2. The definition of the pyredcom functions which can be used using the dataclasses.

References
----------
* Koch, S., Kohrs, F., Lahmann, P., Bissinger, T., Wendschuh, S., Benndorf, D., ... & Klamt, S. (2019).
  RedCom: A strategy for reduced metabolic modeling of complex microbial communities and its
  application for analyzing experimental datasets from anaerobic digestion.
  <i>PLoS computational biology, 15(2)</i>, e1006759.
  [doi:10.1371/journal.pcbi.1006759](https://doi.org/10.1371/journal.pcbi.1006759)
"""
# IMPORT SECTION
# External modules
import cobra
import copy
from dataclasses import dataclass
from typing import Dict, List, Tuple


# DATACLASS DEFINITIONS SECTION
@dataclass
class SingleModel:
    """Dataclass for the description of a single cobra model.

    Description
    ----------
    An instance of this dataclass contains a cobrapy model itself,
    as well as additional information regarding its objective and
    its inputs and outputs.

    Use case
    ----------
    In pyredcom, SingleModel instances are intended to be used
    within Community instances. The additional information of SingleModel
    instances is then used for the generation of the actual community model.

    Member variables
    ----------
    * cobra_model: cobra.Model ~ The cobrapy model
    * species_abbreviation: str ~ A short name for the species that is represented by cobra_model.
      It must not contain underscores!
    * objective_reaction_id: str ~ The reaction ID of the cobra_model's objective (only single-reaction
      objectives are allowed!).
    * exchange_reaction_id_prefix: str ~ The prefix of the cobra_model's exchange reactions. E.g., in
      BIGG models, "EX_" is the prefix for exchange reactions.
    * input_metabolite_ids: List[str] ~ A list of metabolite IDs of metabolites which can be taken up
      by the cobra_model. At least all essential medium metabolite should be listed here.
    * output_metabolite_ids: List[str] ~ A list of metabolite IDs of metabolites which can be taken up
      by the cobra_model. At least all essential medium metabolite should be listed here.
    * model_metabolite_to_exchange_id_mapping: Dict[str, str] ~ A dictionary containing the names of
      all input and output metabolite IDs as keys, and the corresponding metabolite ID names in the
      exchange compartment as children. The children must not already contain the exchange compartment's
      ID itself as this ID will be added at a later stage after its definition in the Community instance.

    Usage example
    ----------
    Suppose we defined a BIGG-compliant core model of <i>E. coli</i>, and
    we set the exchange compartment metabolite names just as the metabolite
    IDs themselves without the addition of the underscore plus compartment name, then
    a SingleModel definition could look like this:
    <pre>
    single_ec_core_model = SingleModel(
        cobra_model=e_coli_core_model,
        objective_reaction_id="BIOMASS_ec_core",
        exchange_reaction_id_prefix="EX_",
        input_metabolite_ids=["glc__D_e", "o2_e", "h_e", "h2o_e"],
        output_metabolite_ids=["ac_e", "h2o_e", "h_e"],
        model_metabolite_to_exchange_id_mapping={
            "glc__D_e": "glc__D",
            "o2_e": "o2",
            "h_e": "h",
            "h2o_e": "h2o",
            "ac_e": "ac",
        }
    )
    </pre>
    """
    cobra_model: cobra.Model
    species_abbreviation: str
    objective_reaction_id: str
    exchange_reaction_id_prefix: str
    input_metabolite_ids: List[str]
    output_metabolite_ids: List[str]
    model_metabolite_to_exchange_id_mapping: Dict[str, str]


@dataclass
class Community:
    """Dataclass for the definition of a community of SingleModel instances.

    Description
    ----------
    An instance of Community contains single SingleModel instances which describe
    the single-organism constituents of the community. In addition, it contains
    information of the community's exchange compartment and what is allowed for
    uptake and intake with the community.

    Use case
    ----------
    In pyredcom, Community instances are directly used in order to create a
    community model according to the RedCom method.

    Member variables
    ----------
    * single_models: List[SingleModel] ~ The list of single models which shall be
      merged into a RedCom community model.
    * exchange_compartment_id: str ~ The name extension of the RedCom-analogous
      exchange compartment metabolites. One common name would be "_exchg".
    * exchange_reaction_id_prefix: str ~ The prefix for the community's (not the
      single organisms's) exhange reactions, i.e. the metabolites that the
      community can have as intake or outtake. One common prefix would be
      "EX_C_" (where the "C_" stands for community :-).
    * input_metabolite_ids: List[str] ~ A list of IDs of exchange compartment
      metabolites which can be taken as input by the community. This list of IDs must not contain
      the exchange_compartment_id.
    * output_metabolite_ids: List[str] ~ A list of IDs of exchange compartment
      metabolites which can be secreted by the community. This list of IDs must not contain
      the exchange_compartment_id.

    Usage example
    ----------
    Suppose we want to create a community of two core models - one from *Escherichia coli* and
    one from *Vibrio natriegens* - and we want to name the exchange compartment "exchg" and
    set the prefix of the community's exchange metabolites to "EX_C" and *E. coli* has the
    input metabolites glc__D, o2, h and h2o, and *V. natriegens* the additional input
    metabolite biotin, and both habe the output metabolites ac and co2,
    then we could set the Community instance as follows:
    <pre>
    community_example = Community(
        single_models = [ecoli_core, vnatriegens_core],
        exchange_compartment_id = "exchg",
        exchange_reaction_id_prefix = "EX_C",
        input_metabolite_ids = ["glc__D", "o2", "h", "h2o", "ac", "biotin"],
        output_metabolite_ids = ["ac", "co2"]
    )
    </pre>
    """
    single_models: List[SingleModel]
    exchange_compartment_id: str
    exchange_reaction_id_prefix: str
    input_metabolite_ids: List[str]
    output_metabolite_ids: List[str]


# FUNCTION DEFINITIONS SECTION
def community_model_fba_summary(model: cobra.Model, exchange_reaction_id_prefix: str = "EX_C_",
                                optimization_title: str = "FBA") -> Tuple[cobra.Solution, List[Dict[str, int]]]:
    """Performs an FBA (or any kind of modified optimization) on a pyredcom-generated community model and prints the results.

    Description
    ----------
    The model-described optimization problem is solved using cobrapy.Model's optimize() function. Then,
    the in and out fluxes of the whole community (i.e., what the community takes and secretes) are printed,
    as well as the species-internal in and out fluxes and the occurring species.


    Return values
    ----------
    Two values are returned:

    1. A cobra.Solution instance containing the FBA solution. This instance can be used as
       any other cobrapy FBA result.
    2. A listed (i.e., put into []) dictionary containing the occurrence or non-occurrence of species in the
       calculated optimization solution. E.g., if a species with the community ID "ecoli" occurs with growth
       in the solution and another community species calles "vnat" does not occur in it, the resulting value
       would be [{"ecoli1": 1, "vnat": 0}]. "1" stands for "occurrence", "0" for "non-occurrence".


    Arguments
    ----------
    * model: cobra.Model ~ The pyredcom-generated RedCom community model as cobrapy Model instance.
    * exchange_reaction_id_prefix: str = "EX_C_" ~ The exchange reaction ID prefix as given during the
      generation of the model's Community instance.
    * optimization_title: str = "FBA"~ Will be printed at the beginning of the FBA summary as
      f"==SUMMARY OF {optimization_title} RESULT=="
    """
    # Optimize for the model's given optimization problem :D
    fba_solution = model.optimize()

    # Get community in and out fluxes text
    community_exchange_reaction_ids = [x.id for x in model.reactions if x.id.startswith(exchange_reaction_id_prefix)]
    community_in_fluxes = "COMMUNITY IN FLUXES:"
    community_out_fluxes = "\nCOMMUNITY OUT FLUXES:"
    for exchange_reaction_id in community_exchange_reaction_ids:
        flux = round(fba_solution.fluxes[exchange_reaction_id], 3)
        description = "\n" + exchange_reaction_id.replace(exchange_reaction_id_prefix, "") + ": " + str(flux)
        if flux < 0:
            community_in_fluxes += description
        elif flux > 0:
            community_out_fluxes += description

    # Get sepecies-internal in and out fluxes text
    species_exchange_reaction_ids = [x.id for x in model.reactions if x.id.startswith("EXCHG_")]
    species_in_fluxes = "\nSPECIES-INTERNAL IN FLUXES:"
    species_out_fluxes = "\nSPECIES-INTERNAL OUT FLUXES:"
    active_organisms: List[str] = []
    inactive_organisms: List[str] = []
    for exchange_reaction_id in species_exchange_reaction_ids:
        flux = round(fba_solution.fluxes[exchange_reaction_id], 3)
        description = "\n" + exchange_reaction_id.replace("EXCHG_", "") + ": " + str(flux)
        if flux < 0:
            species_in_fluxes += description
            active_organisms.append(exchange_reaction_id.split("_")[1])
        elif flux > 0:
            species_out_fluxes += description
            active_organisms.append(exchange_reaction_id.split("_")[1])
        else:
            inactive_organisms.append(exchange_reaction_id.split("_")[1])
    active_organisms = list(set(active_organisms))
    inactive_organisms = list(set(inactive_organisms))

    # Get objective solution text
    objective = f"\nObjective {str(model.objective.expression)} has the value...\n{str(fba_solution.objective_value)}"

    # Create organism activity dictionary and corresponding output text
    organism_occurence_dictionary: Dict[str, int] = {}
    active_organisms_text = f"\nActive organisms: "
    inactive_organisms_text = f"\nInactive organisms: "
    for active_organism in active_organisms:
        organism_occurence_dictionary[active_organism] = 1
        active_organisms_text += "\n* " + active_organism
    for inactive_organism in inactive_organisms:
        if inactive_organism not in organism_occurence_dictionary.keys():
            organism_occurence_dictionary[inactive_organism] = 0
            inactive_organisms_text += "\n* " + inactive_organism

    # Print the optimization results :D
    print(f"===SUMMARY OF {optimization_title}===")
    print(community_in_fluxes)
    print(community_out_fluxes)
    print(species_in_fluxes)
    print(species_out_fluxes)
    print(objective)
    print(active_organisms_text)
    print(inactive_organisms_text)

    return fba_solution, [organism_occurence_dictionary]


def generate_community_cobra_model(community: Community) -> cobra.Model:
    """Creates a RedCom community cobrapy model from a pyredcom Community instance.

    Description
    ----------
    This functions uses all the extra information that is given to the Community
    instance's member variables - and its SingleModel member variables - in order
    to generate a community model with an exchange compartment, species-specific
    exchange reactions as well as community-wide exchange reactions.

    Return value
    ----------
    A community model (according to the RedCom method) in the form of a cobrapy
    Model instance.

    Arguments
    ----------
    * community: Community ~ The Community instance describing the total community.
    """
    # Check that an actual community is given
    if len(community.single_models) <= 1:
        print("ERROR: Less than 2 models in given Community instance!")
        return

    # Check that no underscores are in the organism IDs
    for single_model in community.single_models:
        if "_" in single_model.species_abbreviation:
            print(f"ERROR: Underscore in the given species abbreviation {single_model.species_abbreviation} D:")
            return

    for single_model in community.single_models:
        # Check for exchange reactions with standard name
        exchange_reaction_ids = [x.id for x in single_model.cobra_model.reactions
                                 if x.id.startswith(single_model.exchange_reaction_id_prefix)]
        for exchange_reaction_id in exchange_reaction_ids:
            reaction = single_model.cobra_model.reactions.get_by_id(exchange_reaction_id)
            reaction.lower_bound = 0
            reaction.upper_bound = 0

        # Rename metabolites
        for metabolite in single_model.cobra_model.metabolites:
            metabolite.id += "_" + single_model.species_abbreviation

        # Rename reactions
        for reaction in single_model.cobra_model.reactions:
            reaction.id += "_" + single_model.species_abbreviation

    # Merge models
    merged_cobra_model = copy.deepcopy(community.single_models[0].cobra_model)
    for single_model in community.single_models[1:]:
        merged_cobra_model.merge(single_model.cobra_model, inplace=True)

    # Create exchange compartment
    merged_cobra_model.compartments[community.exchange_compartment_id] = "Exchange compartment"
    for single_model in community.single_models:
        for model_metabolite_id in single_model.model_metabolite_to_exchange_id_mapping.keys():
            community_model_metabolite_id = model_metabolite_id + "_" + single_model.species_abbreviation
            community_model_metabolite = merged_cobra_model.metabolites.get_by_id(community_model_metabolite_id)

            exchange_metabolite_id = single_model.model_metabolite_to_exchange_id_mapping[model_metabolite_id]
            community_model_exchange_metabolite_id = exchange_metabolite_id+"_"+community.exchange_compartment_id
            community_model_metabolite_ids = [x.id for x in merged_cobra_model.metabolites]
            if community_model_exchange_metabolite_id in community_model_metabolite_ids:
                exchange_metabolite = merged_cobra_model.metabolites.get_by_id(community_model_exchange_metabolite_id)
            else:
                exchange_metabolite = cobra.Metabolite(
                    id=exchange_metabolite_id+"_"+community.exchange_compartment_id,
                    compartment=community.exchange_compartment_id
                )

            exchange_reaction_id = "EXCHG_"+single_model.species_abbreviation+"_"+model_metabolite_id+"_to_"+exchange_metabolite_id
            exchange_reaction = cobra.Reaction(id=exchange_reaction_id)
            if model_metabolite_id in single_model.input_metabolite_ids:
                exchange_reaction.lower_bound = -1000
                if model_metabolite_id in single_model.output_metabolite_ids:
                    exchange_reaction.upper_bound = 1000
                else:
                    exchange_reaction.upper_bound = 0
            elif model_metabolite_id in single_model.input_metabolite_ids:
                exchange_reaction.lower_bound = 0
                exchange_reaction.upper_bound = 1000
            exchange_reaction.add_metabolites({
                community_model_metabolite: -1,
                exchange_metabolite: 1
            })
            merged_cobra_model.add_reactions([exchange_reaction])

    # Input metabolites
    for input_metabolite_id in community.input_metabolite_ids:
        community_model_metabolite = merged_cobra_model.metabolites.get_by_id(input_metabolite_id+"_"+community.exchange_compartment_id)
        input_reaction = cobra.Reaction(
            id=community.exchange_reaction_id_prefix+input_metabolite_id+"_"+community.exchange_compartment_id,
            lower_bound=-1000,
            upper_bound=0
        )
        input_reaction.add_metabolites({
            community_model_metabolite: -1,
        })
        merged_cobra_model.add_reactions([input_reaction])

    # Output metabolites
    for output_metabolite_id in community.output_metabolite_ids:
        community_model_metabolite_id = output_metabolite_id + "_" + community.exchange_compartment_id
        community_model_metabolite = merged_cobra_model.metabolites.get_by_id(community_model_metabolite_id)
        output_reaction_id = community.exchange_reaction_id_prefix+output_metabolite_id+"_"+community.exchange_compartment_id
        community_model_reaction_ids = [x.id for x in merged_cobra_model.reactions]
        if output_reaction_id in community_model_reaction_ids:
            reaction = merged_cobra_model.reactions.get_by_id(output_reaction_id)
            reaction.upper_bound = 1000
        else:
            output_reaction = cobra.Reaction(
                id=output_reaction_id,
                lower_bound=0,
                upper_bound=1000
            )
            output_reaction.add_metabolites({
                community_model_metabolite: -1,
            })
            merged_cobra_model.add_reactions([output_reaction])

    # Community biomass metabolite
    community_biomass_metabolite = cobra.Metabolite(
        id="community_biomass",
        compartment=community.exchange_compartment_id
    )

    objective_reaction_ids = [x.objective_reaction_id+"_"+x.species_abbreviation for x in community.single_models]
    for objective_reaction_id in objective_reaction_ids:
        objective_reaction = merged_cobra_model.reactions.get_by_id(objective_reaction_id)
        objective_reaction.add_metabolites({
            community_biomass_metabolite: 1
        })

    # Setting of new objective reaction
    merged_objective_reaction = cobra.Reaction(
        id="COMMUNITY_GROWTH",
        lower_bound=0,
        upper_bound=1000
    )
    merged_objective_reaction.add_metabolites({
        community_biomass_metabolite: -1,
    })
    merged_cobra_model.add_reactions([merged_objective_reaction])
    merged_cobra_model.objective = "COMMUNITY_GROWTH"

    return merged_cobra_model


def minimal_species_redcom_fba(model: cobra.Model, fixed_growth_rate: float,
                               target_reaction_id: str = "", target_reaction_value: float = 0.0,
                               previous_solutions: List[Dict[str, int]] = []) -> Tuple[cobra.Solution, List[Dict[str, int]]]:
    """Performs a minimal species RedCom FBA with the given model and - optionally - the given target reaction flux.

    Description
    ----------
    This function changes the RedCom FBA's objective by searching a solution for which the *minimal* amount of
    species is needed in order to fulfill the given community growth rate and - optionally - the target reaction flux.

    The whole optimization process is called the "minimal species RedCom FBA" and is defined as a MILP using optlang.
    In comparison to the original RedCom FBA formulation (Koch et al., 2019)

    Return values
    ----------
    The same 2 values as this module's community_model_fba_summary(), the following variables and constraints are added:

    1. For each of the community's organism, an integer variable which can be either 0 or 1 is defined. "0" stands for
        "this organism does *not* plays a role in the solution", "1" stands for "this organism plays a role in the solution".
        Let us call each of these variables S_s, where "S" stands for "species indicator", and "s" for the specific species
        it is associated to.
    2. For each S_s, the following new constraint is added: 0 <= S_s - F_s <= 1, where F_s stands for the fraction variable
       of the associated species (see Koch et al., 2019). In other words, if S_s is 0 (i.e., the species does not occur),
       the equation is 0 <= 0 - F_s <= 1. Since F_s is in the range [0;1], the species' fraction F_s must be 0 if S_s
       is also 0 in order to fulfill the newly added constraint.
    3. The objective function was redefined as the *minimization* of the sum of all S_s. In other words, the minimal amount
       of species needed in order to fulfill the given constraints is searched.

    Arguments
    ----------
    * model: cobra.Model ~ The pyredcom-generated community cobrapy model.
    * fixed_growth_rate: float ~ The fixed grwoth rate for all species of the modelled community.
    * target_reaction_id: str = "" ~ A target reaction which can be set to a given value using target_reaction_value.
      Is "" if no target reaction shall be set.
    * target_reaction_value: float = 0.0 ~ If target_reaction_id != "", this is the flux value to which the reaction
      will be set.
    * previous_solutions: List[Dict[str, int]] = [] ~ A list of previous species solutions. All these solutions
      will be disallowed for the current round of the minimal specied RedCom FBA.
    """
    # If given, set the target reaction to the given flux value
    # TODO: Change to yield
    model = copy.deepcopy(model)
    if target_reaction_id != "":
        model.reactions.get_by_id(target_reaction_id).lower_bound = target_reaction_value
        model.reactions.get_by_id(target_reaction_id).upper_bound = target_reaction_value

    # Set fraction and species constraints
    community_biomass_metabolite = model.metabolites.get_by_id("community_biomass")
    organism_id_biomass_reaction_mapping: Dict[str, str] = {}
    for reaction in community_biomass_metabolite.reactions:
        if reaction.id == "COMMUNITY_GROWTH":
            continue

        organism_id = reaction.id.split("_")[-1]
        organism_id_biomass_reaction_mapping[organism_id] = reaction.id

        species_variable =\
            model.problem.Variable(
                name=organism_id+"_integer",
                lb=0,
                ub=1,
                type="integer"
            )
        model.add_cons_vars(species_variable)

        fraction_variable =\
            model.problem.Variable(
                name=organism_id,
                lb=0.0,
                ub=1.0
            )
        model.add_cons_vars(fraction_variable)
        variable_constraint = model.problem.Constraint(
            model.reactions.get_by_id(reaction.id).flux_expression - species_variable*fixed_growth_rate,
            lb=0.0,
            ub=0.0
        )
        model.add_cons_vars(variable_constraint)

        species_constraint = model.problem.Constraint(
            species_variable - fraction_variable,
            lb=0.0,
            ub=1.0
        )
        model.add_cons_vars(species_constraint)

    # Split all organism-specific reactions to irreversible reactions
    # in order to be able to set the fraction variable
    model = split_reversible_organism_reactions(model, organism_id_biomass_reaction_mapping)
    #
    for reaction in model.reactions:
        # In some models, the upper bound is set to "inf". Since this can cause trouble with some
        # solvers, it is changed to 1000 here.
        if reaction.upper_bound > 1000:
            reaction.upper_bound = 1000

        # Retrieve the reaction's organism ID
        organism_id = reaction.id.split("_")[-1]
        if organism_id not in organism_id_biomass_reaction_mapping.keys():
            continue

        # Set lower bound constraint with fraction variable
        original_lower_bound = reaction.lower_bound
        reaction.lower_bound = 0
        fraction_variable = model.variables[organism_id]
        fraction_constraint = model.problem.Constraint(
            model.reactions.get_by_id(reaction.id).flux_expression - original_lower_bound*fraction_variable,
            name=reaction.id+"_lower_bound_with_fraction_variable",
            lb=0,
            ub=2000
        )
        model.add_cons_vars(fraction_constraint)

        # Set upper bound constraint with fraction variable
        original_upper_bound = reaction.upper_bound
        fraction_variable = model.variables[organism_id]
        fraction_constraint = model.problem.Constraint(
            original_upper_bound*fraction_variable - model.reactions.get_by_id(reaction.id).flux_expression,
            name=reaction.id+"_upper_bound_with_fraction_variable",
            lb=0,
            ub=2000
        )
        model.add_cons_vars(fraction_constraint)

    organism_ids = list(organism_id_biomass_reaction_mapping.keys())
    fraction_summation = model.variables[organism_ids[0]]
    for organism_id in organism_ids[1:]:
        fraction_summation += model.variables[organism_id]
    fraction_sum_constraint = model.problem.Constraint(
        -1.0 + fraction_summation,
        name="fraction_sum_equal_one",
        lb=0,
        ub=0
    )
    model.add_cons_vars(fraction_sum_constraint)

    # Set previous result constraints
    for previous_solution in previous_solutions:
        organism_ids = list(previous_solution.keys())
        if len(organism_ids) == 0:
            continue

        first_organism = organism_ids[0]
        variable = model.variables[first_organism+"_integer"]
        integer_cut = previous_solution[first_organism] * variable
        integer_cut += (-1) * variable
        if len(organism_ids) > 1:
            for organism_id in organism_ids[1:]:
                variable = model.variables[organism_id+"_integer"]
                integer_cut += previous_solution[organism_id] * variable
                integer_cut += (-1) * variable

    if len(previous_solutions) > 0:
        fraction_sum_constraint = model.problem.Constraint(
            integer_cut,
            name="integer_cut_previous_solutions",
            ub=-1
        )
        model.add_cons_vars([fraction_sum_constraint])

    # Set growth rate
    model.reactions.get_by_id("COMMUNITY_GROWTH").lower_bound = fixed_growth_rate
    model.reactions.get_by_id("COMMUNITY_GROWTH").upper_bound = fixed_growth_rate

    # Set minimal organism count objective
    organism_variable_summation = model.variables[organism_ids[0]+"_integer"]
    for organism_id in organism_ids[1:]:
        organism_variable_summation += model.variables[organism_id+"_integer"]
    organism_objective = model.problem.Objective(
        organism_variable_summation,
        direction="min"
    )
    model.objective = organism_objective

    # Perform FBA :D
    return community_model_fba_summary(model, optimization_title="REDCOM FBA")


def redcom_fba(model: cobra.Model, fixed_growth_rate: float) -> Tuple[cobra.Solution, List[Dict[str, int]]]:
    """Run a linearized RedCom FBA.

    Description
    ----------
    Performs the RedCom FBA as described in (Koch et al., 2019; see reference section in the module's main
    description).

    Return values
    ----------
    The same 2 values as this module's community_model_fba_summary()

    Arguments
    ----------
    * model: cobra.Model ~ The RedCom-based cobrapy community model.
    * fixed_growth_rate: float ~ The fixed growth rate of the whole community.
    """
    model = copy.deepcopy(model)
    community_biomass_metabolite = model.metabolites.get_by_id("community_biomass")

    # Set fraction constraints
    organism_id_biomass_reaction_mapping: Dict[str, str] = {}
    for reaction in community_biomass_metabolite.reactions:
        if reaction.id == "COMMUNITY_GROWTH":
            continue

        organism_id = reaction.id.split("_")[-1]
        organism_id_biomass_reaction_mapping[organism_id] = reaction.id

        fraction_variable =\
            model.problem.Variable(
                name=organism_id,
                lb=0.0,
                ub=1.0
            )
        model.add_cons_vars(fraction_variable)
        variable_constraint = model.problem.Constraint(
            model.reactions.get_by_id(reaction.id).flux_expression - model.variables[organism_id]*fixed_growth_rate,
            lb=0.0,
            ub=0.0
        )
        model.add_cons_vars(variable_constraint)

    model = split_reversible_organism_reactions(model, organism_id_biomass_reaction_mapping)
    # In some models, the upper bound is set to "inf". Since this can cause trouble with some
    # solvers, it is changed to 1000 here.
    for reaction in model.reactions:
        if reaction.upper_bound > 1000:
            reaction.upper_bound = 1000

        organism_id = reaction.id.split("_")[-1]
        if organism_id not in organism_id_biomass_reaction_mapping.keys():
            continue

        # Lower bound fraction variable constraint
        original_lower_bound = reaction.lower_bound
        reaction.lower_bound = 0
        fraction_variable = model.variables[organism_id]
        fraction_constraint = model.problem.Constraint(
            model.reactions.get_by_id(reaction.id).flux_expression - original_lower_bound*fraction_variable,
            name=reaction.id+"_lower_bound_with_fraction_variable",
            lb=0,
            ub=2000
        )
        model.add_cons_vars(fraction_constraint)

        # Upper bound fraction variable constraint
        original_upper_bound = reaction.upper_bound
        fraction_variable = model.variables[organism_id]
        fraction_constraint = model.problem.Constraint(
            original_upper_bound*fraction_variable - model.reactions.get_by_id(reaction.id).flux_expression,
            name=reaction.id+"_upper_bound_with_fraction_variable",
            lb=0,
            ub=2000
        )
        model.add_cons_vars(fraction_constraint)

    # Fraction summation constraint
    organism_ids = list(organism_id_biomass_reaction_mapping.keys())
    fraction_summation = model.variables[organism_ids[0]]
    for organism_id in organism_ids[1:]:
        fraction_summation += model.variables[organism_id]
    fraction_sum_constraint = model.problem.Constraint(
        -1.0 + fraction_summation,
        name="fraction_sum_equal_one",
        lb=0,
        ub=0
    )
    model.add_cons_vars(fraction_sum_constraint)

    # Set growth rate
    model.reactions.get_by_id("COMMUNITY_GROWTH").lower_bound = fixed_growth_rate
    model.reactions.get_by_id("COMMUNITY_GROWTH").upper_bound = fixed_growth_rate

    # Perform FBA :D
    return community_model_fba_summary(model, optimization_title="MINIMAL SPECIES REDCOM FBA")


def split_reversible_organism_reactions(model: cobra.Model, organism_id_biomass_reaction_mapping: Dict[str, str]) -> cobra.Model:
    """Splits all reversible reactions of a RedCom community model which do not contain exchange compartment metabolites.

    Description
    ----------
    This function returns a model in which all organismic reversible reactions (i.e., minimal flux < 0) are split into non-reversible
    "forward" and "reverse" reactions. An "organismic" reaction is one which occurs without interaction with metabolites of the
    community's exchange compartment.
    E.g., an organismic reaction named "CBD_ecoli" with minimal flux -50 and maximal flux 1000 would be
    separated into the reactions

    1. CBD_fwd_ecoli with flux range [0;1000] and
    2. CBD_rev_ecoli with flux range [0;50]

    Return value
    ----------
    The modified cobrapy model as cobra.Model instance.

    Arguments
    ----------
    * model: cobra.Model ~ A pyredcom-generated community cobrapy model.
    * organism_id_biomass_reaction_mapping: Dict[str, str] ~ A dictionary with the community model's organism IDs
      as keys, and the corresponding biomass reaction IDs as values
    """
    for reaction in model.reactions:
        # Skip irreversible reactions
        if reaction.lower_bound >= 0:
            continue
        # Skip biomass reactions
        if reaction.id in organism_id_biomass_reaction_mapping.values():
            continue

        # Get biomass ID of the reaction using the pyredcom standard naming scheme
        organism_id = reaction.id.split("_")[-1]
        # Skip reaction if it is not part of the organisms (i.e., it is either an ignored
        # organism or part of the exchange compartment)
        if organism_id not in organism_id_biomass_reaction_mapping.keys():
            continue

        # Create new reaction as described in main comment
        new_reaction = copy.deepcopy(reaction)
        original_lower_bound = reaction.lower_bound
        reaction_id_split = reaction.id.split("_")
        reaction.id = "_".join(reaction_id_split[:-1]) + "_fwd_" + reaction_id_split[-1]
        new_reaction.id = "_".join(reaction_id_split[:-1]) + "_rev_" + reaction_id_split[-1]
        reaction.lower_bound = 0
        new_reaction.lower_bound = 0
        new_reaction.upper_bound = -original_lower_bound

        # Reverse direction of products and educts in reverse reaction
        new_reaction_metabolites_copy = copy.deepcopy(new_reaction.metabolites)
        for key in list(new_reaction_metabolites_copy.keys()):
            new_reaction_metabolites_copy[key] *= -2

        # Add new reaction to model
        new_reaction.add_metabolites(new_reaction_metabolites_copy)
        model.add_reactions([new_reaction])

    return model
