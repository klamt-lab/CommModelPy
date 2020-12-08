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
commodelpy. commodelpy's source code is intersected in two parts:

1. The definition of all dataclasses which describe single models and a community.
2. The definition of the commodelpy functions which can be used using the dataclasses.

The definitions of the community structure are inspired by the ones used in the RedCom
publication (Koch et al., 2019).

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
from typing import Dict, List


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
    In commodelpy, SingleModel instances are intended to be used
    within Community instances. The additional information of SingleModel
    instances is then used for the generation of the actual community model.

    Member variables
    ----------
    * cobra_model: cobra.Model ~ The cobrapy model
    * species_abbreviation: str ~ A short name for the species that is represented by cobra_model.
      It must not contain underscores!
    * objective_reaction_id: str ~ The reaction ID of the cobra_model's objective (only single-reaction
      objectives are allowed!).
    * exchange_reaction_id_prefix: str ~ The prefix of the original cobra_model's exchange reactions. E.g., in
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
    Suppose we defined a BiGG-compliant core model of <i>E. coli</i>, and
    we set the exchange compartment metabolite names just as the metabolite
    IDs themselves without the addition of the underscore plus compartment name, then
    a SingleModel definition could look like this:
    <pre>
    single_ec_core_model = SingleModel(
        cobra_model=e_coli_core_model,
        species_abbreviation="ecoli1",
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
    In commodelpy, Community instances are directly used in order to create a
    community model according to the method described in commodelpy's publication.

    Member variables
    ----------
    * single_models: List[SingleModel] ~ The list of single models which shall be
      merged into a community model.
    * exchange_compartment_id: str ~ The name extension of the
      exchange compartment metabolites. One common name would be "_exchg" *which is necessary
      for the usage with the ASTHERISC package*.
    * exchange_reaction_id_prefix: str ~ The prefix for the community's (not the
      single organisms's) exchange reactions, i.e. the metabolites that the
      community can have as intake or outtake. One common prefix would be
      "EX_C_" (where the "C_" stands for community :-) *which is necessary
      for the usage with the ASTHERISC package*.
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


"""
# FUNCTION DEFINITIONS SECTION
def community_model_fba_summary(model: cobra.Model, exchange_reaction_id_prefix: str = "EX_C_",
                                optimization_title: str = "FBA") -> Tuple[cobra.Solution, List[Dict[str, int]]]:
    Performs an FBA (or any kind of modified optimization) on a commodelpy-generated community model and prints the results.

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
    * model: cobra.Model ~ The commodelpy-generated community model as cobrapy Model instance.
    * exchange_reaction_id_prefix: str = "EX_C_" ~ The exchange reaction ID prefix as given during the
      generation of the model's Community instance.
    * optimization_title: str = "FBA"~ Will be printed at the beginning of the FBA summary as
      f"==SUMMARY OF {optimization_title} RESULT=="

    # Optimize for the model's given optimization problem :D
    fba_solution = model.optimize()

    # Get community in and out fluxes text
    community_exchange_reaction_ids = [
        x.id for x in model.reactions if x.id.startswith(exchange_reaction_id_prefix)]
    community_in_fluxes = "COMMUNITY IN FLUXES:"
    community_out_fluxes = "\nCOMMUNITY OUT FLUXES:"
    for exchange_reaction_id in community_exchange_reaction_ids:
        flux = round(fba_solution.fluxes[exchange_reaction_id], 3)
        description = "\n" + \
            exchange_reaction_id.replace(
                exchange_reaction_id_prefix, "") + ": " + str(flux)
        if flux < 0:
            community_in_fluxes += description
        elif flux > 0:
            community_out_fluxes += description

    # Get sepecies-internal in and out fluxes text
    species_exchange_reaction_ids = [
        x.id for x in model.reactions if x.id.startswith("EXCHG_")]
    species_in_fluxes = "\nSPECIES-INTERNAL IN FLUXES:"
    species_out_fluxes = "\nSPECIES-INTERNAL OUT FLUXES:"
    active_organisms: List[str] = []
    inactive_organisms: List[str] = []
    for exchange_reaction_id in species_exchange_reaction_ids:
        flux = round(fba_solution.fluxes[exchange_reaction_id], 3)
        description = "\n" + \
            exchange_reaction_id.replace("EXCHG_", "") + ": " + str(flux)
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
"""


def generate_community_model_with_no_growth(community: Community, fractions: Dict[str, float], biomass_reactions: Dict[str, str] = {}) -> cobra.Model:
    """Creates a cobrapy-compatible community model from a commodelpy Community instance.

    Description
    ----------
    This function uses all information that is given to the Community
    instance's member variables - and its SingleModel member variables - in order
    to generate a community model with an exchange compartment, species-specific
    exchange reactions as well as community-wide exchange reactions as defined
    by the input and output metabolites of the respective Community instance
    and its SingleSpecies instances.

    Herein, the model-specific fractions of all single organisms can be given and for every species-associated
    function with a minimal/maximal flux which is not equal to inf/-inf/0, the
    minimal/maximal flux (denoted as v) is set to v*fraction_of_the_reaction's_species.
    Since only the growth rate is now free, this linearizes the 'problem' of solving
    a balanced-growth community model for growth, and an FBA on this model will show
    the theoretical optimal growth with the given fixed organism fractions.

    Return value
    ----------
    A community model with fixed species ratios and no growth in the form of a cobrapy
    Model instance. This form of the community model - with its fixed organism
    ratios - can be e.g. used with the ASTHERISC package.

    Arguments
    ----------
    * community: Community ~ The Community instance describing the total community. *In order to run with the ASTHERISC package,
      the community's exchange compartment ID must be set to "exchg" and the community's exchange reaction ID prefix to "EX_C"*.
    * fractions: Dict[str, float]: float ~ A dictionary containing all species names of
      the community as keys, and their fractions of the total community biomass as
      values. For reasonable calculations, the sum of all fractions should be 1.0.
    * biomass_reactions: Dict[str, str] ~ If not {}, it must denote the biomass reaction IDs (values) for each
      of the community's species (keys). For these reactions, the mock biomass metabolite community_biomass_$SPECIES_NAME
      will be added as product and the objective function community_biomass_reaction will be able to consume it.
      If it is {}, no objective function will be created and the mock biomass metabolites will not be added to the species biomass reaction.
      Is {} by default.
    """
    # Get number of SingleModel instances in Community instance
    num_single_models = len(community.single_models)
    # Check that an actual community is given
    if num_single_models <= 1:
        raise ValueError(
            "ERROR: Less than 2 models in given Community instance!")
    # Check that the right amount of fractions is given
    if len(fractions.keys()) != num_single_models:
        raise ValueError(
            f"ERROR: Number of given fractions ({len(fractions)}) does not match number of single models ({num_single_models})!")

    # Go through each SingleModel and change their cobra models
    for single_model in community.single_models:
        species_fraction = fractions[single_model.species_abbreviation]

        # Check that no underscores are in the organism IDs
        if "_" in single_model.species_abbreviation:
            raise ValueError(
                f"ERROR: Underscore in the given species abbreviation {single_model.species_abbreviation} D:")

        # Rename metabolites
        for metabolite in single_model.cobra_model.metabolites:
            metabolite.id += "_" + single_model.species_abbreviation

        # Rename reactions and change bounds
        for reaction in single_model.cobra_model.reactions:
            reaction.id += "_" + single_model.species_abbreviation

            if reaction.upper_bound != float("inf"):
                reaction.upper_bound *= species_fraction
            if reaction.lower_bound != -float("inf"):
                reaction.lower_bound *= species_fraction
            if (reaction.upper_bound == float("inf")) and (species_fraction == 0):
                reaction.upper_bound = 0
            if (reaction.lower_bound == -float("inf")) and (species_fraction == 0):
                reaction.lower_bound = 0

        # Delete standard exchange reactions with default exchange reaction ID prefix
        standard_exchanges = [x for x in single_model.cobra_model.reactions
                              if x.id.startswith(single_model.exchange_reaction_id_prefix)]
        single_model.cobra_model.remove_reactions(standard_exchanges)

    # Merge single models into a huge one
    merged_model = copy.deepcopy(community.single_models[0].cobra_model)
    for single_model in community.single_models[1:]:
        merged_model.merge(single_model.cobra_model, inplace=True)

    # Add whole community exchanges
    exchange_metabolite_ids = list(
        set(community.input_metabolite_ids + community.output_metabolite_ids))
    for exchange_metabolite_id in exchange_metabolite_ids:
        # Set reaction instance
        reaction = cobra.Reaction(id=community.exchange_reaction_id_prefix+exchange_metabolite_id+"_"+community.exchange_compartment_id,
                                  name="Community exchange for "+exchange_metabolite_id)

        # Set reaction bounds
        is_input = exchange_metabolite_id in community.input_metabolite_ids
        if is_input:
            reaction.lower_bound = -float("inf")
        else:
            reaction.lower_bound = 0
        is_output = exchange_metabolite_id in community.output_metabolite_ids
        if is_output:
            reaction.upper_bound = float("inf")
        else:
            reaction.upper_bound = 0

        # Add metabolite to reaction
        exchange_metabolite = cobra.Metabolite(exchange_metabolite_id+"_"+community.exchange_compartment_id,
                                               name="Exchange compartment metabolite "+exchange_metabolite_id,
                                               compartment="exchange")
        reaction.add_metabolites({
            exchange_metabolite: -1,
        })

        # Add reaction to model
        merged_model.add_reactions([reaction])

    # Add mock community biomass metabolites for the ASTHERISC package's species recognition
    biomass_metabolite = cobra.Metabolite(id="community_biomass",
                                          name="Mock community biomass metabolite for ASTHERISC package species recognition",
                                          compartment="exchg")
    merged_model.add_metabolites([biomass_metabolite])

    # Add single species <-> exchange compartment exchanges
    for single_model in community.single_models:
        exchange_metabolite_ids = list(
            set(single_model.input_metabolite_ids + single_model.output_metabolite_ids))
        exchange_metabolite_ids = [
            x+"_"+single_model.species_abbreviation for x in exchange_metabolite_ids]

        # Add mock community biomass metabolites for the ASTHERISC package's species recognition
        biomass_metabolite = cobra.Metabolite(id="community_biomass_"+single_model.species_abbreviation,
                                              name="Mock community biomass metabolite for ASTHERISC package species recognition",
                                              compartment="exchg")
        merged_model.add_metabolites([biomass_metabolite])

        for exchange_metabolite_id in exchange_metabolite_ids:
            exchange_compartment_metabolite_id = single_model.model_metabolite_to_exchange_id_mapping[exchange_metabolite_id.replace(
                "_"+single_model.species_abbreviation, "")]

            # Set reaction instance
            reaction = cobra.Reaction(id="EXCHG_"+single_model.species_abbreviation+"_"+exchange_metabolite_id.replace("_"+single_model.species_abbreviation, "")+"_to_"+exchange_compartment_metabolite_id,
                                      name=f"Exchange for {exchange_metabolite_id} from single species {single_model.species_abbreviation} to exchange compartment")

            # Set reaction bounds
            # is_input = exchange_compartment_metabolite_id in single_model.input_metabolite_ids
            is_input = exchange_metabolite_id.replace("_"+single_model.species_abbreviation, "") in single_model.input_metabolite_ids
            if is_input:
                reaction.lower_bound = -float("inf")
            else:
                reaction.lower_bound = 0
            # is_output = exchange_compartment_metabolite_id in single_model.output_metabolite_ids
            is_output = exchange_metabolite_id.replace("_"+single_model.species_abbreviation, "") in single_model.output_metabolite_ids
            if is_output:
                reaction.upper_bound = float("inf")
            else:
                reaction.upper_bound = 0

            # Add metabolites to reaction
            internal_metabolite = merged_model.metabolites.get_by_id(
                exchange_metabolite_id)
            exchange_compartment_metabolite = merged_model.metabolites.get_by_id(
                exchange_compartment_metabolite_id+"_"+community.exchange_compartment_id)
            reaction.add_metabolites({
                internal_metabolite: -1,
                exchange_compartment_metabolite: 1,
            })

            # Add reaction to model
            merged_model.add_reactions([reaction])

    return merged_model


def create_community_model_with_balanced_growth(community: Community, growth_rate: float) -> cobra.Model:
    """Creates a combined community model with an stoichiometric-matrix-integrated balanced growth approach.

    Description
    ----------
    *Aim of this function*

    This method is aimed to generate a community model (as described in the community argument of
    this function) with a fixed growth rate, where this growth rate is the same for each of the community's
    species, in a way in which this growth rate constraint is directly integrated into the community model
    itself so that this model can be used e.g. with common Flux Balace Analysis methods and funtions.

    *Basic background behind this function*

    In community models, both the fraction of single species on the total community
    biomass as well as the growth rate of each single species can be optimized, which would lead to a computaionally
    quite complex bilinear optimization. Since balanced growth (i.e., no species grows faster or slower than other species
    so that no species is outcompeted in the long term) requires a fixed growth rate, only the optimization of the species
    fractions remains and, therefore, a computationally much more efficient linear optimiztaion.

    As described e.g. in (Koch et al., 2019), the resulting minimal flux of an irreversible reaction j in species i with a minimal flux greater than 0
    is min_flux_j*f_i where f is the species's fraction, and, consequently, the resulting maximal flux (if the maximum is smalelr than inf) is max_flux_j*f_i. Usually,
    a special Flux Balance Analysis function introducing the fraction variables as added variables has to be written in order
    to run FBA with balanced growth.

    This unconvenient introduction of extra fraction variables can be ommitted by looking at the meaning of a species's fraction:
    Since the fractions of single species (if we give them the index i) f_i = µ_i/µ_community, where µ is the growth rate,
    and µ_community is fixed the only thing one to look at is µ_i, i.e. the "biomass contribution" of each species which
    is equal to their fraction.

    Hence, the simplified approach used heirin, which implicitly integrates the fraction variables into the stoichiometric
    matrix, is as follows:

    * For reactions with a minimal flux > 0: A new pseudo-metabolite is integrated as where it is produced by the reactions itself with the stoichiometry 1,
      and this metabolite is consumed by the reaction's species biomass reaction with the stoichiometry min_flux/fixed_growth_rate. Additionally, a
      pseudo-reaction consuming the pseudo-metabolite is added, too.
    * For reactions with a maximal flux < inf: A new pseudo-metabolite is integrated as where it is produced by the reactions itself with the stoichiometry 1,
      and this metabolite is consumed by the reaction's species biomass reaction with the stoichiometry max_flux/fixed_growth_rate. Additionally, a
      pseudo-reaction producing the pseudo-metabolite is added, too.


    *Principle of balanced growth constraint integration used herein*

    1. All metabolites and reactions of the single species (as defined in the SingleModel instances of the given Community instance)
       are renamed with the addition of "_"+species_abbreviation (species_abbreviation as defined in the SingleModel instance) at
       the end of their IDs.
    2. A new pseudo-metabolite called "COMMUNITY_BIOMASS" is added as product with stoichiometry 1 to all biomass reactions
       of the single species.
    3. All single species biomass reactions get the lower flux bound of 0 and upper flux bound of infinite.
    4. A new pseudo-reaction called "COMMUNITY_BIOMASS" is generated which has the upper and lower (i.e., fixed) flux equal
       to the given growth rate. This pseudo-reaction has the pseudo-metabolite "COMMUNITY_BIOMASS" as educt with a
       stoichiometry of 1.
    5. The pseudo-reactions (called "Rsnake_UPPER_" for maximal flux bounds and "Rsnake_LOWER_" for minimal flux bounds)
       and pseudo-metabolites (called "Msnake_UPPER_" for maximal flux bounds and "Msnake_LOWER_" for minimal flux bounds)
       are introduced into the model as described in the previous section.
    6. Exchange reactions between the single species comparments and the exchange compartment for all input/output
       metabolites from the respective SingleSpecies instances, as well as from the exchange
       compartment to the environment for all input/output metabolites of the Community instance. All previous standard
       exchanges of the single models (usually starting with "EX_") are deleted.

    The principle used herein also means that if another growth rate is to be tested, a new model has to be generated for this
    specific growth rate.

    Return value
    ----------
    A community model in the form of a cobrapy
    Model instance. This form of the community model - with its fixed organism
    ratios - can be e.g. used with the ASTHERISC package.

    Arguments
    ----------
    * community: Community ~
    * growth_rate: float ~
    """
    # Get number of SingleModel instances in Community instance
    num_single_models = len(community.single_models)
    # Check that an actual community is given
    if num_single_models <= 1:
        raise ValueError(
            "ERROR: Less than 2 models in given Community instance!")

    # Dictionary for later biomass introduction
    organism_id_biomass_reaction_id_mapping: Dict[str, str] = {}

    # Go through each SingleModel and change their cobra models
    for single_model in community.single_models:
        # Check that no underscores are in the organism IDs
        if "_" in single_model.species_abbreviation:
            raise ValueError(
                f"ERROR: Underscore in the given species abbreviation {single_model.species_abbreviation} D:")

        # Rename metabolites
        for metabolite in single_model.cobra_model.metabolites:
            metabolite.id += "_" + single_model.species_abbreviation

        # Rename reactions and change bounds
        for reaction in single_model.cobra_model.reactions:
            reaction.id += "_" + single_model.species_abbreviation

        # Delete standard exchange reactions with default exchange reaction ID prefix
        standard_exchanges = [x for x in single_model.cobra_model.reactions
                              if x.id.startswith(single_model.exchange_reaction_id_prefix)]
        single_model.cobra_model.remove_reactions(standard_exchanges)

        # Add biomass reaction to dictionary
        organism_id_biomass_reaction_id_mapping[single_model.species_abbreviation] = single_model.objective_reaction_id + \
            "_" + single_model.species_abbreviation
        # Standardize single-organism biomass reaction bounds
        biomass_reaction = single_model.cobra_model.reactions.get_by_id(
            organism_id_biomass_reaction_id_mapping[single_model.species_abbreviation])
        biomass_reaction.lower_bound = 0
        biomass_reaction.upper_bound = float("inf")

    # Merge single models into a huge one
    merged_model = copy.deepcopy(community.single_models[0].cobra_model)
    for single_model in community.single_models[1:]:
        merged_model.merge(single_model.cobra_model, inplace=True)

    # Add community biomass metabolite
    community_biomass_metabolite = cobra.Metabolite(
        id="COMMUNITY_BIOMASS", name="Community biomass metabolite", compartment="exchg")
    merged_model.add_metabolites([community_biomass_metabolite])
    for biomass_reaction_id in organism_id_biomass_reaction_id_mapping.values():
        biomass_reaction = merged_model.reactions.get_by_id(
            biomass_reaction_id)
        biomass_reaction.add_metabolites({
            community_biomass_metabolite: 1,
        })

    # Add whole community exchanges
    exchange_metabolite_ids = list(
        set(community.input_metabolite_ids + community.output_metabolite_ids))
    for exchange_metabolite_id in exchange_metabolite_ids:
        # Set reaction instance
        reaction = cobra.Reaction(id=community.exchange_reaction_id_prefix+exchange_metabolite_id+"_"+community.exchange_compartment_id,
                                  name="Community exchange for "+exchange_metabolite_id)

        # Set reaction bounds
        is_input = exchange_metabolite_id in community.input_metabolite_ids
        if is_input:
            reaction.lower_bound = -float("inf")
        else:
            reaction.lower_bound = 0
        is_output = exchange_metabolite_id in community.output_metabolite_ids
        if is_output:
            reaction.upper_bound = float("inf")
        else:
            reaction.upper_bound = 0

        # Add metabolite to reaction
        exchange_metabolite = cobra.Metabolite(exchange_metabolite_id+"_"+community.exchange_compartment_id,
                                               name="Exchange compartment metabolite "+exchange_metabolite_id,
                                               compartment="exchange")
        reaction.add_metabolites({
            exchange_metabolite: -1,
        })

        # Add reaction to model
        merged_model.add_reactions([reaction])

    # Split reversible reactions
    merged_model = split_reversible_organism_reactions(
        merged_model, organism_id_biomass_reaction_id_mapping)

    # Add community biomass reaction and set it to the given growth rate
    community_biomass_reaction = cobra.Reaction(id="COMMUNITY_BIOMASS",
                                                name="Biomass reaction for the whole community")
    community_biomass_reaction.lower_bound = growth_rate
    community_biomass_reaction.upper_bound = growth_rate
    community_biomass_reaction.add_metabolites({
        community_biomass_metabolite: -1,
    })

    # Add minimal and maximal bound constraint
    reaction_ids = [x.id for x in merged_model.reactions]
    for reaction_id in reaction_ids:
        reaction = merged_model.reactions.get_by_id(reaction_id)

        # Check organism ID
        reaction_organism_id = reaction.id.split("_")[-1]
        if reaction_organism_id not in organism_id_biomass_reaction_id_mapping.keys():
            continue

        organism_biomass_reaction = merged_model.reactions.get_by_id(
            organism_id_biomass_reaction_id_mapping[reaction_organism_id])
        # Add maximal bound constraint
        if reaction.upper_bound != float("inf"):
            # Add ~M
            new_metabolite = cobra.Metabolite(id="Msnake_UPPER_"+reaction.id,
                                              name="Upper bound enforcing metabolite for "+reaction.id,
                                              compartment="exchg")
            # Add ~r
            new_reaction = cobra.Reaction(id="Rsnake_UPPER_"+reaction.id,
                                          name="Delivery reaction of upper bound enforcing metabolite for "+reaction.id)
            # Add ~M to original reaction
            reaction.add_metabolites({
                new_metabolite: 1,
            })
            organism_biomass_reaction.add_metabolites({
                new_metabolite: -reaction.upper_bound/growth_rate,
            })
            new_reaction.add_metabolites({
                new_metabolite: 1,
            })
            merged_model.add_reactions([new_reaction])

        # Add minimal bound constraint
        if (reaction.lower_bound != -float("inf")) and (reaction.lower_bound != 0.0):
            # Add ~M
            new_metabolite = cobra.Metabolite(id="Msnake_LOWER_"+reaction.id,
                                              name="Lower bound enforcing metabolite for "+reaction.id,
                                              compartment="exchg")
            # Add ~r
            new_reaction = cobra.Reaction(id="Rsnake_LOWER_"+reaction.id,
                                          name="Delivery reaction of lower bound enforcing metabolite for "+reaction.id)
            # Add ~M to original reaction
            reaction.add_metabolites({
                new_metabolite: 1,
            })
            organism_biomass_reaction.add_metabolites({
                new_metabolite: -reaction.lower_bound/growth_rate,
            })
            new_reaction.add_metabolites({
                new_metabolite: -1,
            })
            merged_model.add_reactions([new_reaction])

    # Add single species <-> exchange compartment exchanges
    for single_model in community.single_models:
        exchange_metabolite_ids = list(
            set(single_model.input_metabolite_ids + single_model.output_metabolite_ids))
        exchange_metabolite_ids = [
            x+"_"+single_model.species_abbreviation for x in exchange_metabolite_ids]

        # Add mock community biomass metabolites for the ASTHERISC package's species recognition
        biomass_metabolite = cobra.Metabolite(id="community_biomass_"+single_model.species_abbreviation,
                                              name="Mock community biomass metabolite for ASTHERISC package species recognition",
                                              compartment="exchg")
        merged_model.add_metabolites([biomass_metabolite])

        for exchange_metabolite_id in exchange_metabolite_ids:
            exchange_compartment_metabolite_id = single_model.model_metabolite_to_exchange_id_mapping[exchange_metabolite_id.replace(
                "_"+single_model.species_abbreviation, "")]

            # Set reaction instance
            reaction = cobra.Reaction(id="EXCHG_"+single_model.species_abbreviation+"_"+exchange_metabolite_id.replace("_"+single_model.species_abbreviation, "")+"_to_"+exchange_compartment_metabolite_id,
                                      name=f"Exchange for {exchange_metabolite_id} from single species {single_model.species_abbreviation} to exchange compartment")

            # Set reaction bounds
            is_input = exchange_compartment_metabolite_id in community.input_metabolite_ids
            if is_input:
                reaction.lower_bound = -float("inf")
            else:
                reaction.lower_bound = 0
            is_output = exchange_compartment_metabolite_id in community.output_metabolite_ids
            if is_output:
                reaction.upper_bound = float("inf")
            else:
                reaction.upper_bound = 0

            # Add metabolites to reaction
            internal_metabolite = merged_model.metabolites.get_by_id(
                exchange_metabolite_id)
            exchange_compartment_metabolite = merged_model.metabolites.get_by_id(
                exchange_compartment_metabolite_id+"_"+community.exchange_compartment_id)
            reaction.add_metabolites({
                internal_metabolite: -1,
                exchange_compartment_metabolite: 1,
            })

            # Add reaction to model
            merged_model.add_reactions([reaction])

    # Set merged model's objective to community biomass
    merged_model.objective = "COMMUNITY_BIOMASS"

    return merged_model


def split_reversible_organism_reactions(model: cobra.Model, organism_id_biomass_reaction_mapping: Dict[str, str]) -> cobra.Model:
    """Splits all reversible reactions of a community model which do not contain exchange compartment metabolites.

    Description
    ----------
    This function returns a model in which all organismic reversible reactions (i.e., minimal flux < 0) are split into non-reversible
    "forward" and "reverse" reactions. An "organismic" reaction is one which occurs without interaction with metabolites of the
    community's exchange compartment.
    E.g., an organismic reaction named "ABC_ecoli" with minimal flux -50 and maximal flux 1000 would be
    separated into the reactions

    1. ABC_forward_ecoli with flux range [0;1000] and
    2. ABC_reverse_ecoli with flux range [0;50]

    Return value
    ----------
    The modified cobrapy model as cobra.Model instance.

    Arguments
    ----------
    * model: cobra.Model ~ A commodelpy-generated community cobrapy model.
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

        # Get biomass ID of the reaction using the commodelpy standard naming scheme
        organism_id = reaction.id.split("_")[-1]
        # Skip reaction if it is not part of the organisms (i.e., it is either an ignored
        # organism or part of the exchange compartment)
        if organism_id not in organism_id_biomass_reaction_mapping.keys():
            continue

        # Create new reaction as described in main comment
        new_reaction = copy.deepcopy(reaction)
        original_lower_bound = reaction.lower_bound
        reaction_id_split = reaction.id.split("_")
        reaction.id = "_".join(
            reaction_id_split[:-1]) + "_forward_" + reaction_id_split[-1]
        new_reaction.id = "_".join(
            reaction_id_split[:-1]) + "_reverse_" + reaction_id_split[-1]
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
