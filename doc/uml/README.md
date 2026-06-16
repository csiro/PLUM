# PLUM UML Diagrams

This directory contains UML diagrams documenting the PLUM architecture and execution flow.

## Available Diagrams

### Behavioral Diagrams (Execution Flow)

#### 1. reftest_sequence.puml

**UML Sequence Diagram: Reference Test Execution Flow**

This comprehensive sequence diagram traces the complete execution path from the bash script `reftest.sh` through the PLUM executable down to C++ class-level interactions. It emphasizes **temporal ordering** and **execution flow**. It shows:

**Bash Script Layer:**
- Command-line argument parsing
- Test file reading and processing
- Command execution via `eval`
- Output validation with `diff`

**Main Executable Layer (plum.cpp):**
- Configuration and option parsing (Opts, Config)
- Parameter setup (Params)
- Random seed initialization (Rand)
- Scenario loading (metabolites, reactions, experiments)
- Solver selection logic (CTS, INT, INCR, DUMMY)
- LP/MILP implementation instantiation (Gurobi, HiGHS, lp_solve)

**C++ Class Architecture:**
- **Scenario**: Metabolic network container with metabolites, reactions, compartments
- **Reaction**: Stoichiometric transformations with flux bounds and objective costs
- **Metabolite**: Compounds with supply/demand rates across compartments
- **GapSolver** (abstract): Base solver interface
  - LPSolver: Continuous LP formulation (flux minimization)
  - IntSolver: Integer programming (reaction selection)
  - IncrSolver: Incremental solver (progressive addition)
- **LPSolverImp** (abstract): LP/MILP implementation interface
  - GrbLPSolverImp: Gurobi backend
  - HighsLPSolverImp: HiGHS backend
  - LpsLPSolverImp: lp_solve backend
- **Solution**: Result container with flux values and statistics

**Key Interactions Shown:**
- LP/MILP problem formulation (stoichiometry, bounds, objectives)
- Solver optimization loop
- Solution extraction and statistics calculation
- Multiple output formats (flux, model, metabolite balance)

**Best for:** Understanding the chronological flow of execution, debugging, tracing method calls

---

#### 2. reftest_communication.puml

**UML Communication Diagram: Object Relationships and Message Flow**

This communication diagram shows the same system as the sequence diagram but emphasizes **object relationships** and **structural connections** rather than temporal ordering. It shows:

**Key Features:**
- All objects involved in the test execution
- Communication links between objects
- Numbered messages indicating sequence (1-86)
- Bidirectional relationships and data flow
- Multiple output responsibilities per object

**Structural Emphasis:**
- Which objects communicate with which
- Message aggregation on links
- Object roles and responsibilities
- Package/layer organization implied

**Differences from Sequence Diagram:**
- Spatial layout shows relationships
- Messages numbered but not temporally ordered
- Multiple messages shown on single link
- Easier to see "who talks to whom"
- Better for understanding architecture

**Best for:** Understanding system architecture, identifying coupling, seeing communication patterns, design reviews

---

#### 3. reftest_communication_compact.puml

**UML Communication Diagram (Compact): Key Structural Relationships**

A simplified communication diagram that focuses on the most important structural relationships and message paths. It uses **package notation** to group related components.

**Key Features:**
- Package-based organization (7 logical layers)
- Inheritance relationships shown with <|--
- Aggregated message descriptions
- Design pattern annotations
- Communication flow summary

**Packages:**
1. **Test Layer**: User, reftest.sh
2. **Main Executable**: plum main()
3. **Configuration**: Opts, Config, Params, Rand
4. **Data Model**: Scenario, Metabolite, Reaction, Experiment
5. **Solver Framework**: GapSolver (interface), LPSolver, IntSolver, IncrSolver
6. **LP/MILP Implementation**: LPSolverImp (interface), GrbLPSolverImp, HighsLPSolverImp, LpsLPSolverImp
7. **Results**: Solution, MultiSol

**Design Patterns Identified:**
- Strategy Pattern (solver selection, backend selection)
- Facade Pattern (PlumMain orchestration)
- Composite Pattern (Scenario aggregation)

**Best for:** High-level architecture overview, identifying design patterns, presentations, documentation

---

### Structural Diagrams (Class Architecture)

#### 4. class_solver_hierarchy.puml

**UML Class Diagram: Solver Framework and LP/MILP Implementations**

This comprehensive class diagram shows the complete solver architecture with inheritance hierarchies and implementation strategies. It shows:

**Solver Framework (GapSolver Hierarchy):**
- Abstract GapSolver base class
- LPSolver: Continuous LP formulation (flux minimization)
- IntSolver: Integer programming (MILP with binary use variables)
- IncrSolver: Incremental solver (cost-based reaction addition)
- DummySolver: Testing solver
- LnsMxSolver, MathHeur: Advanced solvers

**LP/MILP Implementation Layer (LPSolverImp Hierarchy):**
- Abstract LPSolverImp interface
- GrbLPSolverImp: Gurobi backend (commercial)
- HighsLPSolverImp: HiGHS backend (open source)
- LpsLPSolverImp: lp_solve backend (open source)

**Key Features:**
- Complete method signatures with parameters and return types
- Public/protected/private access modifiers
- Detailed member variables and their purposes
- Inheritance relationships (<|-- notation)
- Composition relationships (*-- notation)
- Design pattern annotations (Strategy, Composite)
- Algorithm notes (LP formulation, MILP formulation, incremental algorithm)

**Best for:** Understanding solver architecture, implementing new solvers, refactoring solver code, algorithm comparison

---

#### 5. class_data_model.puml

**UML Class Diagram: Metabolic Network Data Model**

This class diagram details the core data structures representing metabolic networks for flux balance analysis. It shows:

**Core Classes:**
- Scenario: Top-level container (metabolites, reactions, experiments)
- Metabolite: Biochemical compounds with supply/demand rates
- Reaction: Stoichiometric transformations with flux bounds
- Experiment: Growth conditions (media composition)
- SupplyResidual: Supply and residual rates per metabolite

**Solution Classes:**
- Solution: Flux distribution with objective values
- MultiSol: Multiple experiment solutions aggregator

**Configuration:**
- Params: Solver parameters and configuration

**Utilities:**
- PathFinder: Reachability analysis
- DualVals: Dual values from LP solution

**Key Features:**
- Complete attribute lists with types
- Comprehensive method signatures
- Relationship cardinalities (1-to-many, etc.)
- File I/O methods (PLD format)
- Compartment model annotations
- Special metabolite/reaction types
- Stoichiometry representation

**Best for:** Understanding data structures, implementing file parsers, adding new features to data model, database schema design

---

#### 6. class_system_overview.puml

**UML Class Diagram: Complete System Architecture**

This high-level overview diagram shows all major components and their interactions across the entire PLUM system. It shows:

**Package Organization:**
1. **lime (Utility Library)**: Config, Opts, Rand, TimeKeeper, Dig, Displayable
2. **Data Model**: Scenario, Metabolite, Reaction, Experiment, Params
3. **Solver Framework**: GapSolver hierarchy (LPSolver, IntSolver, IncrSolver)
4. **LP/MILP Backends**: LPSolverImp hierarchy (Gurobi, HiGHS, lp_solve)
5. **Solution**: Solution, MultiSol
6. **Utilities**: PathFinder, DualVals, enums (StatusEnum, Flavour)

**Key Features:**
- Package-level organization
- Simplified class interfaces (key methods only)
- All major relationships shown
- System data flow legend
- Design pattern identification
- Command-line usage example

**Best for:** System overview, onboarding, presentations, high-level architecture documentation

---

## Diagram Comparison

| Aspect | Sequence | Communication | Comm. Compact | Class Diagrams |
|--------|----------|---------------|----------------|
| **Emphasis** | Time/order | Relationships | Structure | Static structure |
| **Layout** | Temporal | Spatial | Hierarchical |
| **Messages** | Ordered arrows | Numbered on links | Aggregated |
| **Best for** | Execution flow | Object interactions | Architecture |
| **Detail level** | High | Medium | Low | High |
| **Complexity** | High | Medium | Low | Medium-High |
| **Use case** | Debugging | Design review | Overview | Architecture |
| **File size** | ~12KB | ~7KB | ~5KB | ~10-15KB |

**When to use which:**
- **Sequence**: Debugging, understanding execution, following code flow
- **Communication**: Architecture review, identifying dependencies, refactoring
- **Communication Compact**: Presentations, onboarding, high-level documentation
- **Class Diagrams**: Understanding structure, implementing features, refactoring classes

## Generating Diagrams

### Prerequisites

Install PlantUML and its dependencies:

```bash
# Ubuntu/Debian
sudo apt-get install plantuml

# Or download PlantUML JAR
wget https://github.com/plantuml/plantuml/releases/download/v1.2024.7/plantuml-1.2024.7.jar

# Requires Java runtime
sudo apt-get install default-jre

# For better rendering, install Graphviz
sudo apt-get install graphviz
```

### Generate PNG Images

Using system PlantUML:
```bash
cd doc/uml

# Generate specific diagrams
plantuml reftest_sequence.puml
plantuml reftest_communication.puml
plantuml reftest_communication_compact.puml
plantuml class_solver_hierarchy.puml
plantuml class_data_model.puml
plantuml class_system_overview.puml
```

Or generate all at once:
```bash
cd doc/uml
plantuml *.puml
```

Using PlantUML JAR:
```bash
cd doc/uml
java -jar plantuml-1.2024.7.jar *.puml
```

This generates:
- `reftest_sequence.png` (sequence diagram)
- `reftest_communication.png` (detailed communication diagram)
- `reftest_communication_compact.png` (compact communication diagram)
- `class_solver_hierarchy.png` (solver framework class diagram)
- `class_data_model.png` (data model class diagram)
- `class_system_overview.png` (system overview class diagram)

### Generate SVG (Scalable Vector Graphics)

```bash
plantuml -tsvg *.puml
```

This generates SVG files for all diagrams (recommended for documentation)

### Generate Other Formats

```bash
# PDF
plantuml -tpdf reftest_sequence.puml

# LaTeX
plantuml -tlatex reftest_sequence.puml

# ASCII art (text-based, no graphics needed)
plantuml -ttxt reftest_sequence.puml

# EPS (for publications)
plantuml -teps reftest_sequence.puml
```

### Batch Processing

Generate all diagrams at once:
```bash
plantuml doc/uml/*.puml
```

### Real-time Preview

Some editors provide PlantUML preview:
- **VS Code**: Install "PlantUML" extension
- **IntelliJ IDEA**: Built-in PlantUML support
- **Vim**: plantuml-syntax plugin
- **Emacs**: plantuml-mode

### Online Rendering

For quick viewing without installation:
1. Copy the contents of any `.puml` file
2. Visit: https://www.plantuml.com/plantuml/uml/
3. Paste and view the rendered diagram

**Recommended for each diagram type:**
- **Sequence**: Best viewed as PNG or in IDE with large vertical space
- **Communication**: Best viewed as SVG for pan/zoom
- **Communication Compact**: Best viewed as PNG or in presentations
- **Class Diagrams**: Best viewed as SVG for clarity, or PNG for documentation

## Understanding the Diagrams

### Reading the Sequence Diagram

The sequence diagram flows **top-to-bottom**, **left-to-right**:

1. **User invokes reftest.sh** (top)
2. **Bash script processes test file** (initialization)
3. **Loop through test commands** (main test loop)
4. **Execute plum command** (transition to C++)
5. **C++ initialization** (config, params, scenario)
6. **Solver setup** (flavour and type selection)
7. **Problem formulation** (LP/MILP setup)
8. **Optimization** (solver execution)
9. **Results extraction** (solution statistics)
10. **Output generation** (files written)
11. **Return to bash** (validation with diff)
12. **Test completion** (bottom)

### Diagram Elements

- **Boxes at top**: Participants (scripts, executables, classes)
- **Vertical lines**: Lifelines (participant activity duration)
- **Arrows**: Messages/method calls between participants
- **Activation bars**: When a participant is actively processing
- **Notes**: Additional context and explanations
- **Alt/Loop boxes**: Control flow (conditionals, loops)

### Key Sections

- **Test Initialization**: Argument parsing, file validation
- **Test Execution Loop**: Reading and executing commands
- **PLUM Executable Invocation**: Transition from bash to C++
- **Scenario Loading**: PLD file parsing, network construction
- **Solver Setup**: Implementation and type selection
- **Problem Solving**: LP/MILP formulation and optimization
- **Results Processing**: Statistics and solution analysis
- **Output Generation**: Writing various result formats
- **Test Validation**: Comparing output to reference files

## Use Cases

### For Developers

- **Understanding code flow**: Follow the execution path from command to result
- **Debugging**: Identify where data transforms occur
- **Integration**: See how components interact
- **Onboarding**: Quick overview of system architecture

### For Documentation

- **Technical specifications**: Formal system behavior documentation
- **Publications**: High-quality diagrams for papers
- **Presentations**: Visual aid for explaining PLUM architecture
- **Design reviews**: Discuss system design with stakeholders

### For Testing

- **Test coverage**: Identify all code paths exercised by tests
- **Integration points**: See where components connect
- **Validation flow**: Understand how results are verified

### Reading the Communication Diagrams

Communication diagrams emphasize **object relationships** and use **spatial layout**:

**Key Elements:**
- **Objects**: Shown as rectangles with names and roles
- **Links**: Lines between objects showing communication paths
- **Messages**: Numbered arrows on links (sequence from reftest_sequence.puml)
- **Packages**: Grouped components (compact version)

**Reading Strategy:**

1. **Identify layers** (in compact version):
   - Test Layer (top)
   - Main Executable
   - Configuration
   - Data Model
   - Solver Framework
   - LP/MILP Implementation
   - Results (bottom)

2. **Follow message numbers** (1-86):
   - Same numbering as sequence diagram
   - Shows execution order within structural layout

3. **Observe relationships**:
   - Thick links = heavy communication
   - Bidirectional arrows = mutual dependency
   - <|-- notation = inheritance/implementation

4. **Note patterns**:
   - Star topology = facade/orchestrator
   - Chain = data pipeline
   - Multiple inheritance = strategy pattern

**Communication vs Sequence:**
- **Sequence**: "What happens when?"
- **Communication**: "Who talks to whom?"
- Use together for complete understanding

## Diagram Maintenance

When modifying PLUM code, update the UML diagrams to reflect:
- New classes or participants
- Changed method signatures
- Additional control flow paths
- New solver types or implementations
- Modified file formats or outputs

Keep diagrams synchronized with code to maintain documentation accuracy.

## PlantUML Syntax Quick Reference

```plantuml
' Comments start with apostrophe
participant "Name" as Alias    ' Define participant
Alias -> Other: message()      ' Synchronous call
Alias --> Other: return value  ' Return/response
activate Alias                 ' Start activation
deactivate Alias              ' End activation

alt condition                  ' Conditional (if/else)
    ...
else alternative
    ...
end

loop description               ' Loop construct
    ...
end

note left/right of Alias      ' Add notes
    Note text
end note
```

For complete PlantUML syntax: https://plantuml.com/sequence-diagram

## File Organization

```
doc/uml/
├── README.md                              # This file
├── SUMMARY.md                             # Creation summary and feasibility
├── reftest_sequence.puml                  # Sequence diagram (temporal)
├── reftest_communication.puml             # Communication diagram (detailed)
├── reftest_communication_compact.puml     # Communication diagram (compact)
├── reftest_sequence.png                   # Generated images (after plantuml)
├── reftest_communication.png              # ...
├── reftest_communication_compact.png      # ...
└── *.svg                                  # SVG versions (after plantuml -tsvg)
```

## Additional Resources

- **PlantUML Official**: https://plantuml.com/
- **Sequence Diagram Guide**: https://plantuml.com/sequence-diagram
- **UML Specification**: https://www.omg.org/spec/UML/
- **PLUM Documentation**: See `/home/kubeflow/PLUM/CLAUDE.md`
- **Architecture Overview**: See `/home/kubeflow/PLUM/doc/overview/overview.md`

## Contributing Diagrams

When adding new UML diagrams:

1. **Name consistently**: `<component>_<type>.puml` (e.g., `solver_class.puml`)
2. **Include title**: Use `title` directive for diagram heading
3. **Add theme**: Use `!theme plain` for consistent styling
4. **Document**: Add section to this README explaining the diagram
5. **Generate**: Produce PNG/SVG for version control
6. **Validate**: Ensure diagram renders correctly before committing

### Recommended Diagram Types

Already implemented:
- ✅ **Sequence diagrams**: Execution flows (reftest_sequence.puml)
- ✅ **Communication diagrams**: Object relationships (reftest_communication*.puml)
- ✅ **Class diagrams**: Static structure (class_*.puml)

Future additions:
- **Object diagrams**: Instance-level examples
- **Component diagrams**: High-level architecture (executables, libraries)
- **Activity diagrams**: Algorithm flows (reachability analysis, incremental solver)
- **State diagrams**: Solver state machines (if applicable)
- **Deployment diagrams**: System deployment and dependencies

## Questions or Issues

For questions about:
- **PlantUML syntax**: See https://plantuml.com/
- **PLUM architecture**: See `/home/kubeflow/PLUM/CLAUDE.md`
- **Diagram accuracy**: Check source code in `src/` and `include/mosh/`
- **Rendering problems**: Ensure Graphviz is installed
