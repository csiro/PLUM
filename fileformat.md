# Internal File Format

## File has 3 sections

- Compartments (COMPARTMENT)
- Metabolites (MET)
- Reactions (REACTION)

They must appear in that order



### Compartments

```julia
COMPART <id>
```

- id: The compartment ID

- The order of compartments is used to specify values in metabolites and reactions.

  

### Metabolites

```julia
MET <id> <supply-1> <supply-2>... <residual-1> <residual-2> ...
```

- 'id' is the metabolite name
- Supply values specify the rate at which the metabolite is taken from the growth medium
- Residual values specify the rate at which the metabolite is exported from the cell.
- There must be as many supply values as there are compartments, in the order that compartments were specified
- There must be as many residual values as there are compartments, in order

### Reactions

```julia
REACTION <id> <obj-coeff-1> <obj-coeff-2> ... <compound-spec-1> <compound-spec-2> ...
```

- 'id' is the reaction name

- 'obj coeff' gives an objective coefficient for each compartment

- There must be one obj coeff for each compartment, in the order compartments were specified

- 'compound-specs' are as follows

  - `<compound-id> <compartment> <rate>`

  - compound-id is the compound involved
  - 'compartment' is an integer, giving the relative compartment. A transport reaction will have different compartments. Compartment 0 is where reaction products reside, and compartments 1, 2, etc are the source of compounds. With Compartments given in the order OCPE, a reaction with compartment 1 and compartment 0 will be transporting E->P, P->C and/or C->O

## EOF