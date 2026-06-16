# PLUM Class Diagram Guide
## Understanding the Structural Architecture

This guide explains the three class diagrams created for PLUM's structural architecture and how to use them effectively.

---

## Quick Reference

| Diagram | Focus | Level | Best For |
|---------|-------|-------|----------|
| **Solver Hierarchy** | Algorithms & backends | Detailed | Implementation, extending solvers |
| **Data Model** | Network representation | Detailed | Data structures, file I/O |
| **System Overview** | Complete architecture | High-level | Onboarding, presentations |

---

## 1. Solver Hierarchy (class_solver_hierarchy.puml)

### Purpose
Shows the **complete solver framework** with both algorithm strategies and LP/MILP backend implementations.

### What It Contains

**GapSolver Hierarchy (Algorithm Layer):**
```
GapSolver (abstract)
├── LPSolver (continuous LP)
├── IntSolver (integer MILP)
├── IncrSolver (incremental)
├── DummySolver (testing)
├── LnsMxSolver (large neighborhood search)
└── MathHeur (math-heuristic)
```

**LPSolverImp Hierarchy (Backend Layer):**
```
LPSolverImp (abstract interface)
├── GrbLPSolverImp (Gurobi)
├── HighsLPSolverImp (HiGHS)
└── LpsLPSolverImp (lp_solve)
```

### Key Features

**Detailed Method Signatures:**
```cpp
// From LPSolver class
+ solve() : SolutionPtr
+ formulate(exp : Experiment*)
+ do_solve() : SolutionPtr
+ set_target_flux(flux)
+ summary() : string
```

**Member Variables with Types:**
```cpp
// From IntSolver class
- which_formulation_ : int
- imp_ : LPSolverImpPtr
- cts_sol_ : MultiSolPtr
- biomass_obj_mult_ : double
- rounding_iters_ : int
- mip_gap_ : double
```

**Design Patterns Annotated:**
- **Strategy Pattern** (GapSolver): Multiple algorithms
- **Strategy Pattern** (LPSolverImp): Multiple backends
- **Composite Pattern** (IncrSolver): Contains LPSolver

**Algorithm Notes:**
- LP formulation with constraints
- MILP formulation with binary variables
- Incremental algorithm flow

### When to Use

✅ **Implementing a new solver algorithm**
- See what methods to override
- Understand base class responsibilities
- Check member variables needed

✅ **Adding a new LP/MILP backend**
- See LPSolverImp interface contract
- Compare existing implementations
- Understand variable/constraint management

✅ **Understanding solver selection**
- See how -v flag maps to classes
- Understand -V flag backend selection
- See composition relationships

✅ **Debugging solver behavior**
- Trace method calls through hierarchy
- Understand state management
- See interaction between layers

### Example Usage

**Question**: "How do I add support for CPLEX?"

**Answer using diagram**:
1. Find LPSolverImp interface
2. See all methods marked `{abstract}`
3. Note GrbLPSolverImp implementation as example
4. Create CplexLPSolverImp implementing same interface
5. Implement: init_cts(), init_int(), make_flux_var(), add_met_constraint(), optimize(), etc.

---

## 2. Data Model (class_data_model.puml)

### Purpose
Shows the **core data structures** representing metabolic networks for flux balance analysis.

### What It Contains

**Core Classes:**
- **Scenario**: Top-level container
- **Metabolite**: Compounds with properties
- **Reaction**: Transformations with stoichiometry
- **Experiment**: Growth conditions
- **SupplyResidual**: Supply/demand rates

**Solution Classes:**
- **Solution**: Flux distribution
- **MultiSol**: Multiple experiments

**Supporting Classes:**
- **Params**: Configuration parameters
- **PathFinder**: Reachability analysis
- **DualVals**: LP dual values

### Key Features

**Complete Attribute Lists:**
```cpp
// From Metabolite class
- name_ : string
- full_name_ : string
- index_ : size_t
- is_source_ : bool
- is_cycle_met_ : bool
- is_c_source_ : bool
- is_dummy_ : bool
```

**Relationship Cardinalities:**
```
Scenario "1" *-- "many" Metabolite : contains
Scenario "1" *-- "many" Reaction : contains
Reaction "1" o-- "many" Metabolite : uses/makes
```

**File I/O Methods:**
```cpp
+ read_data(filename)
+ read_supply_demand(filename, params)
+ read_react_cost(filename, params, policy)
+ write_to(out)
```

**Domain Annotations:**
- Compartment model (E, P, C, O)
- Special metabolite types
- Reaction classification
- Stoichiometry representation

### When to Use

✅ **Understanding data structures**
- See all attributes and their types
- Understand object relationships
- Find accessor methods

✅ **Implementing file parsers**
- See what data needs to be loaded
- Understand field mappings
- Check validation requirements

✅ **Adding new features to data model**
- See where new attributes fit
- Understand existing patterns
- Check impact on related classes

✅ **Database schema design**
- Map classes to tables
- Identify foreign keys
- Understand normalization needs

### Example Usage

**Question**: "How is stoichiometry stored?"

**Answer using diagram**:
1. Find Reaction class
2. See `mets_ : vector<Metabolite*>` (participating metabolites)
3. See `met_coeff_ : map<Metabolite*, double>` (stoichiometric coefficients)
4. Note methods: `met_coeff(met)`, `uses(met)`, `makes(met)`
5. Understand: coefficient > 0 = product, < 0 = reactant

---

## 3. System Overview (class_system_overview.puml)

### Purpose
Shows the **complete system architecture** with all major components and package organization.

### What It Contains

**Six Major Packages:**

1. **lime (Utility Library)**
   - Config, Opts, Rand, TimeKeeper, Dig, Displayable

2. **Data Model**
   - Scenario, Metabolite, Reaction, Experiment, Params

3. **Solver Framework**
   - GapSolver hierarchy (LPSolver, IntSolver, IncrSolver, DummySolver)

4. **LP/MILP Backends**
   - LPSolverImp hierarchy (Gurobi, HiGHS, lp_solve)

5. **Solution**
   - Solution, MultiSol

6. **Utilities**
   - PathFinder, DualVals, StatusEnum, Flavour

### Key Features

**Simplified Interfaces:**
- Only key methods shown (5-10 per class)
- Focuses on public API
- Emphasizes relationships over details

**Package Organization:**
- Clear separation of concerns
- Logical grouping
- Dependency direction

**System Flow:**
```
Input (PLD files) 
  → Scenario (Data Model)
  → GapSolver (Solver Framework)
  → LPSolverImp (LP/MILP Backend)
  → Solution
  → Output (files)
```

**Design Pattern Summary:**
- Strategy: GapSolver, LPSolverImp
- Factory: Solution creation
- Composite: Scenario aggregation
- Facade: GapSolver interface

### When to Use

✅ **System overview**
- Quick understanding of architecture
- See all major components
- Identify subsystems

✅ **Onboarding new developers**
- First diagram to show
- Understand system organization
- See package dependencies

✅ **Presentations and documentation**
- High-level architecture slide
- Technical reports
- Stakeholder communication

✅ **Planning major changes**
- Identify affected packages
- Understand system boundaries
- Evaluate architectural impact

### Example Usage

**Question**: "What are the major subsystems?"

**Answer using diagram**:
1. See 6 packages clearly delineated
2. Trace data flow: Input → Data Model → Solver → Backend → Solution → Output
3. Identify extension points: GapSolver, LPSolverImp
4. Note utility layer: lime library for common operations

---

## Reading Class Diagrams

### UML Notation Guide

**Access Modifiers:**
```
+ public    : Accessible from anywhere
# protected : Accessible from subclasses
- private   : Accessible only within class
```

**Relationship Types:**
```
<|--  : Inheritance (is-a)
        Example: IntSolver <|-- GapSolver (IntSolver IS-A GapSolver)

*--   : Composition (strong ownership, part-of)
        Example: Scenario *-- Metabolite (Scenario OWNS Metabolites)

o--   : Aggregation (weak ownership, has-a)
        Example: GapSolver o-- Scenario (GapSolver HAS-A Scenario reference)

..>   : Dependency (uses)
        Example: GapSolver ..> Solution (GapSolver CREATES Solution)
```

**Method Notation:**
```
{abstract} method() : ReturnType
           ↑        ↑     ↑
       abstract  method  return type
                  name

+ method(param : Type) : ReturnType
  ↑       ↑       ↑         ↑
public parameter type   return type
```

**Multiplicity:**
```
"1" : Exactly one
"*" : Zero or more
"many" : Multiple
"0..1" : Zero or one
"1..*" : One or more
```

### Reading Strategy

**Top-Down Approach (System Overview → Details):**

1. **Start with System Overview**
   - Get oriented: packages and flow
   - Identify major subsystems
   - Understand high-level relationships

2. **Dive into relevant detail diagram**
   - Solver Hierarchy: For algorithm/backend details
   - Data Model: For structure/storage details

3. **Cross-reference with code**
   - Find class in `include/mosh/`
   - Match methods and members
   - Verify relationships

**Bottom-Up Approach (Class → System):**

1. **Start with specific class**
   - Find in Data Model or Solver Hierarchy
   - Read all attributes and methods
   - Understand responsibilities

2. **Follow relationships**
   - See what it depends on (o--, ..>)
   - See what depends on it (reverse arrows)
   - Understand coupling

3. **Zoom out to System Overview**
   - See how class fits in architecture
   - Identify package membership
   - Understand system context

---

## Practical Workflows

### Implementing a New Feature

**Scenario**: Add support for reaction activation energy

**Workflow**:

1. **Data Model Analysis**
   - Open `class_data_model.puml`
   - Find `Reaction` class
   - Check existing attributes
   - Identify where to add `activation_energy_`

2. **Check Dependencies**
   - See what uses Reaction (Solution, Scenario, GapSolver)
   - Note file I/O methods (`read_data()`)
   - Understand stoichiometry access

3. **Implementation Plan**
   - Add `activation_energy_` to Reaction
   - Add accessor: `activation_energy() : double`
   - Update `read_data()` to parse new field
   - Update `write_to()` to output new field

4. **Verify with System Overview**
   - Check `class_system_overview.puml`
   - Ensure no other subsystems affected
   - Confirm Scenario → Reaction relationship still valid

### Debugging a Solver Issue

**Scenario**: IntSolver returning wrong objective value

**Workflow**:

1. **Solver Hierarchy Investigation**
   - Open `class_solver_hierarchy.puml`
   - Find `IntSolver` class
   - Check `solve()` method (where it's defined)
   - Note member: `mip_gap_ : double`

2. **Trace Implementation**
   - See `IntSolver` uses `imp_ : LPSolverImpPtr`
   - Check which backend in use (GrbLPSolverImp?)
   - Note abstract methods: `get_objective()`, `make_sol()`

3. **Check Formulation**
   - Read algorithm notes in diagram
   - Understand MILP formulation
   - Check biomass_obj_mult_ calculation

4. **Cross-reference with Solution**
   - Open `class_data_model.puml`
   - Find `Solution` class
   - Check `obj_value()`, `rel_obj_value()`, `abs_obj_value()`
   - Understand which objective is being used

### Refactoring

**Scenario**: Extract common LP formulation logic

**Workflow**:

1. **Identify Duplication**
   - Open `class_solver_hierarchy.puml`
   - Compare `LPSolver` and `IntSolver`
   - Note both call `imp_->make_flux_var()`, `imp_->add_met_constraint()`

2. **Design Extraction**
   - Create new class: `FormulationBuilder`
   - Move common methods from LPSolver/IntSolver
   - Update relationships in diagram

3. **Update Diagrams**
   - Add `FormulationBuilder` to Solver Hierarchy
   - Show composition: `LPSolver *-- FormulationBuilder`
   - Update System Overview with new package structure

4. **Validate Architecture**
   - Check System Overview for consistency
   - Ensure package boundaries maintained
   - Verify no circular dependencies introduced

---

## Comparison with Behavioral Diagrams

### Class Diagrams vs Sequence Diagrams

| Aspect | Class Diagram | Sequence Diagram |
|--------|---------------|------------------|
| **Shows** | Static structure | Dynamic behavior |
| **Focus** | Classes, attributes, methods | Messages, order, timing |
| **Time** | No time dimension | Temporal flow |
| **Use for** | Architecture, refactoring | Debugging, tracing |
| **Answers** | "What exists?" | "What happens?" |

**Example**:
- **Class**: "Reaction has a method `met_coeff(met : Metabolite*) : double`"
- **Sequence**: "LPSolver calls `react->met_coeff(met)` during formulation"

### When to Use Both

**Complete Understanding** = Class Diagram + Sequence Diagram

1. **Class Diagram** tells you:
   - What classes exist
   - What methods they have
   - How they're related

2. **Sequence Diagram** tells you:
   - When methods are called
   - In what order
   - With what parameters

**Example workflow**:
1. Read class diagram to understand `IntSolver` structure
2. Read sequence diagram to see how `IntSolver.solve()` executes
3. Return to class diagram to understand helper methods called

---

## Advanced Tips

### Finding Specific Information

**"Where is X stored?"**
→ Data Model diagram, look for X as attribute

**"Who creates Y?"**
→ Look for `..>` dependency arrows pointing to Y

**"What's the inheritance hierarchy of Z?"**
→ Follow `<|--` arrows from Z upward

**"Which classes use W?"**
→ Look for associations (o--, *--, ..>) pointing to W

### Understanding Design Decisions

**Why is there a GapSolver base class?**
- See multiple derived classes (LPSolver, IntSolver, IncrSolver)
- Note pattern: Strategy Pattern for algorithm selection
- Enables runtime selection via command-line flag (-v)

**Why separate LPSolverImp from LPSolver?**
- See LPSolver contains LPSolverImp
- Multiple LPSolverImp implementations (Gurobi, HiGHS, lp_solve)
- Pattern: Strategy Pattern for backend selection
- Separates algorithm (LPSolver) from optimization engine (LPSolverImp)

**Why does Scenario own Metabolite/Reaction but GapSolver references them?**
- Scenario uses composition (*--): Strong ownership, lifecycle management
- GapSolver uses aggregation (o--): Weak reference, doesn't own
- Reason: Scenario is data container, GapSolver is algorithm using that data

### Identifying Extension Points

**To add a new solver algorithm:**
- Extend GapSolver
- Implement `solve() : SolutionPtr`
- Use existing LPSolverImp for backend

**To add a new optimization backend:**
- Extend LPSolverImp
- Implement all abstract methods
- Use with existing GapSolver algorithms

**To add new metabolite properties:**
- Add attribute to Metabolite class
- Add accessor methods
- Update Scenario.read_data() for parsing
- Check Solution output methods if needed

---

## Tools and Rendering

### Recommended Viewing

**Solver Hierarchy:**
- Best as SVG (complex, benefits from zoom)
- Landscape orientation
- Large display recommended

**Data Model:**
- Best as PNG or SVG
- Portrait or landscape
- Good for printing

**System Overview:**
- Best as PNG for presentations
- Fits on single slide/page
- Clear at various sizes

### PlantUML Tips for Class Diagrams

**Hide specific methods:**
```plantuml
class MyClass {
  + important_method()
  - {method} private_method()
}
hide MyClass method private_method
```

**Group related classes:**
```plantuml
together {
  class ClassA
  class ClassB
}
```

**Add custom notes:**
```plantuml
note left of ClassName
  Important information
  about this class
end note
```

---

## Maintenance Guidelines

### Keeping Diagrams Updated

**When to update class diagrams:**

| Change | Update Solver Hierarchy | Update Data Model | Update System Overview |
|--------|------------------------|-------------------|----------------------|
| Add class | ✅ Yes | ✅ Yes | ✅ Yes |
| Add method | ✅ Yes (if public) | ✅ Yes (if key) | ❌ No |
| Add attribute | ✅ Yes | ✅ Yes | ❌ No |
| Change relationship | ✅ Yes | ✅ Yes | ✅ Yes |
| Refactor | ✅ Yes | ✅ Yes | ✅ Yes |
| Bug fix (no structure change) | ❌ No | ❌ No | ❌ No |

### Validation Checklist

Before committing updated diagrams:

- [ ] All class names match source code
- [ ] Inheritance arrows point in correct direction (<|-- = extends)
- [ ] Composition shows ownership (*-- = owns)
- [ ] Aggregation shows reference (o-- = references)
- [ ] Public/protected/private symbols correct (+/#/-)
- [ ] Abstract methods marked with {abstract}
- [ ] Method signatures match implementation
- [ ] Attribute types are accurate
- [ ] Diagram renders without errors
- [ ] Legend/notes are still accurate

---

## Conclusion

The three class diagrams provide **complete structural documentation** for PLUM:

1. **Solver Hierarchy**: Deep dive into algorithms and backends
2. **Data Model**: Complete data structure reference
3. **System Overview**: High-level architecture map

Use them together with behavioral diagrams (sequence, communication) for **complete system understanding**.

**Quick decision tree:**
```
Need to understand structure?
├─ YES → Class diagrams
│   ├─ Whole system? → System Overview
│   ├─ Solvers/algorithms? → Solver Hierarchy
│   └─ Data structures? → Data Model
└─ Need to understand behavior? → Sequence/Communication diagrams
```

---

**See also**:
- `README.md`: Generation instructions
- `DIAGRAM_GUIDE.md`: Sequence vs Communication comparison
- `INDEX.md`: Complete diagram index
- Source code: `include/mosh/*.h`

---

**Created**: 2026-06-15  
**Author**: Claude Code (Anthropic)  
**Project**: PLUM Metabolic Gap-Filling Solver
