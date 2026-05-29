
# Gapfill Project

- Context is finding which reactions take place in a metabolic process 
- a.k.a. Flux Balance Analysis - finding the flux in reactions that balances known inputs with known outputs, or maximises the biomass reaction
- Some reactions are known (or suspected) - we want to find which others are used
- Relies on steady-state reaction
- A database of potential reactions is required

# I/O

Input data

- The rate of use of known input metabolites (e.g. those available in the growth media)
- The rate of production of known output metabolites (usually biomass products)
- A list of possible reactions, with
  - Stoichiometric coefficients giving ratio of use/production of each metabolite
  - The source/destination compartment (Cytosol, Periplasm, External) for each metabolite
  - Notional "cost" of the reaction, reflecting the likelihood that this reaction is used.

Output

- The list of the reactions that are used, along with the flux for each.

# Modelling

GlobalFit (2016) is/was one of the leading method for gap filling. It uses Mixed Integer Linear Programming (MILP)

- Decision variables: 
  - $x_i \in \{0, 1\}, x_i = 1$ if I use reaction $i$ from the database, 0 otherwise.
  - $y_i$ the flux through reaction $i$
- Cost $c_i$ of using reaction $i$ depends on likelihood: small likelihood $\rightarrow$ large cost.
- Data is $a_{ij}$, the stoichiometric coefficient of metabolite $j$ in reaction $i$
- Objective is to minimise the **cost of used reactions** (regardless of the flux)

Model is
$$
\min \sum_i c_i x_i \\
\mbox{subject to} \sum_i a_{ij} y_i = 0 \;\;\; \forall \mbox{ metabolites } j \\
0 \le y_i \le \mbox{UB}_i \; x_i \;\;\; \forall \mbox{ reactions } i \\
x_i \in \{ 0, 1\}
$$


Unfortunately, MILP can be very slow to solve - hours to weeks to arrive at a good solution, potentially years to solve exactly.

# Modelling II

An alternative objective is to minimise the **total flux** used (regardless of the number of reactions)
$$
\min \sum_i c_i y_i \\
\mbox{subject to} \sum_i a_{ij} y_i = 0 \;\;\; \forall \mbox{ metabolites } j \\
0 \le y_i \le {\mbox{UB}}_i \;\;\; \forall \mbox{ reactions } i
$$
This is an non-integer Linear Program (LP). LP can be solved in polynomial time, and in practice is solved in seconds.

Unfortunately we are not the first to think of this. I have just found a paper from 2014 that suggest this approach

- M. Latendresse,  Efficiently gap-filling reaction networks,  *BMC Bioinformatics*, vol. 15, no. 1, p. 225, Jun. 2014, doi: [10.1186/1471-2105-15-225](https://doi.org/10.1186/1471-2105-15-225).)

It's a bit wierd to me that this has not filtered through - a paper from 2021 still uses a MILP technique to solve Community gap-filling problem (2 organisms sharing metabolites)

- D. Giannari, C. H. Ho, and R. Mahadevan,  A gap-filling algorithm for prediction of metabolic interactions in microbial communities,  *PLOS Computational Biology*, vol. 17, no. 11, p. e1009060, Nov. 2021, doi: [10.1371/journal.pcbi.1009060](https://doi.org/10.1371/journal.pcbi.1009060).

# Reachability

- Idea of reachability is a kind of depth analysis. All metabolites available in the growth media are called level 0 (or *seed*) metabolites. Level 1 reactions are those that only use level 0 metabolites as reactants. Level 1 reactions produce level 1 metabolites. At level 2 are all those reactions that can proceed with level 1 or level 0 metabolites, and they produce level 2 metabolites. So level $n$ reactions use only metabolites produced at levels strictly less than $n$.
- Level analysis could help discriminate between possible solutions to a gap-filling scenario - only allowed to choose a reaction if it is "reachable" using already-selected reactions.
- Unfortunately this also turns out to be not a new idea. The concept is called *graph-based producibility* in 
  - C. Frioux,  Investigating host-microbiota cooperation with gap-filling optimization problems, Doctorial Thesis, Université Rennes  2018.
- and used in the computational tool *Meneco*
  - Prigent, S., Frioux, C., Dittami, S. M., Thiele, S., Larhlimi, A., Collet, G., Gutknecht, F., Got, J., Eveillard, D., Bourdon, J., Plewniak, F., Tonon, T., and Siegel, A. (2017). Meneco, a Topology-Based Gap-Filling Tool Applicable to Degraded Genome-Wide Metabolic Networks. PLOS Computational Biology, 13(1):e1005276.
- In a Masters thesis, Thullier (2019) describes a way of incorporating "topological activation" in a hybrid solution framework (although they use a MILP for solution of the FBA )
  - K. Thuillier,  Linear programming for metabolic network completion,  HAL, INRIA, hal-02408003, Dec. 2019.

# Costs

Besides the reactions with their stoichiometric coefficients, the main input data is the "cost" for each equation.

The cost has to reflect the likelihood that the reaction is used. This cost will be very problem-dependent. Considerations in determining the cost include

- There are bench test showing this reaction is used
- There are bench test showing the reaction is a possible match
- I have seen this reaction used in similar situations
- I have seen this reaction used in this species 
- I have seen this reaction used in similar species
- I have seen this reaction used in a vaguely related species
- I have seen this reaction used.

Experimentation will be needed to set the levels for different interpretations of likelihood. 

# Costs II

Just FYI, Latendresse *et al* (2014) use 5 cost levels

1. A cost of 1 for one candidate spontaneous (i.e., nonenzymatic) reaction
2. A cost of 5 for one candidate reaction inside the taxonomic range of the organism
3. A cost of 10 for one candidate reaction of unknown taxonomic range
4. A cost of 15 for one candidate reaction outside the taxonomic range of the organism
5. A gain of 500 for one candidate biomass metabolite 

# First Results

- Results with this MILP and LP models are encouraging
- Can identify main reactions
- However, without Compartments, some reactions are missed. 
  - In particular, modelling the biomass reaction is not accurate


# Compartments

Now investigating use of compartments.

- Reactions occur in Compartments - External, Periplasm, and Cytosol. 

- For our purposes, there are two types of External 
  - External where metabolites are available as *ab initio* reactants
  - External where output products reside. 
- Databases do not distinguish between these externals, but our balance equations have to. 
- We will refer to the output External as Output. We refer to the compartments by there initial (E, P, C, O)

# Example I

A simple reaction in the database is represented like this, for metabolites X and Y

```julia
-1 X(0) +1 Y(0)
```

This means the reaction uses one unit of X to produce one unit of Y in compartment 0. The compartments match, so this is a non-transport reaction. It can happen in any of the compartments - E, P, C or O. We will therefor have to repeat each reaction in each compartment - so the above reaction will lead to the following in our model

```julia
-1 X(E) +1 Y(E)
-1 X(P) +1 Y(P)
-1 X(C) +1 Y(C)
-1 X(O) +1 Y(O)
```

The cost associated with each of these can be different, if there is a convenient way to capture that information.

# Example II

An example transport reaction in the database is represented like this, for metabolites A, B, C

```julia
-1 A(0) -1 B(0) +1 C(1)
```

which means 1 unit of A and 1 unit of B are used to create 1 unit of C, and furthermore, the A and the B are used up in Compartment "0", and C is produced in compartment "1"

This reaction in the forward direction can mean transport from E to P, P to C, or C to O. We will produce 3 reactions for each transport, representing these different cases.

(Just a little detail, we will also have to have a 0-cost transport reaction for every metabolite from E to O. Since E and O are actually the same place, if we have said there is a metabolite in E, we need to allow a reaction in O to use it)

# Example III

- Biomass export will only be modelled in the Output compartment
- All metabolites required for biomass will appear in the biomass reaction, drawing from the Output compartment.

