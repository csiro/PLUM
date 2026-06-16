# PLUM UML Diagram Guide
## Understanding Sequence vs Communication Diagrams

This guide explains the three UML diagrams created for PLUM's reference test execution and when to use each one.

---

## Quick Reference Table

| Diagram | Purpose | Emphasis | Best For | Complexity |
|---------|---------|----------|----------|------------|
| **Sequence** | Execution flow | Temporal ordering | Debugging, code tracing | High |
| **Communication** | Object interactions | Relationships | Design review, refactoring | Medium |
| **Communication Compact** | System structure | Architecture | Overview, presentations | Low |

---

## 1. Sequence Diagram (reftest_sequence.puml)

### Purpose
Shows **when** things happen and in what **order**.

### Visual Characteristics
- **Vertical timeline**: Time flows from top to bottom
- **Horizontal participants**: Objects/classes arranged left to right
- **Lifelines**: Vertical dashed lines showing object lifetime
- **Activation bars**: Rectangles showing when object is processing
- **Arrows**: Messages flowing between participants

### What It Shows Well
✅ **Temporal ordering**: See exact sequence of method calls  
✅ **Nested calls**: See when objects call other objects and return  
✅ **Control flow**: Loops, conditionals, alternatives clearly visible  
✅ **Execution context**: Activation bars show call stack depth  
✅ **Return values**: Dashed arrows show data coming back  

### What It Doesn't Show Well
❌ **Overall structure**: Hard to see "who talks to whom" at a glance  
❌ **Object relationships**: Connections not spatially emphasized  
❌ **Design patterns**: Not obvious from temporal layout  

### Example: Following a Bug

**Scenario**: Biomass flux is zero when it shouldn't be

**How to use the sequence diagram:**
1. Find "biomass_flux()" call (message ~72)
2. Trace backward to see where Solution was created (message ~65)
3. Follow back to solver optimize() call (message ~53)
4. Check formulation messages (messages ~49-52)
5. See what data was passed from Scenario (messages ~42-44)

You can see the **complete call chain** in chronological order.

---

## 2. Communication Diagram (reftest_communication.puml)

### Purpose
Shows **who** talks to **whom** and **how** they interact.

### Visual Characteristics
- **Object boxes**: Components shown as rectangles with roles
- **Communication links**: Lines between objects that communicate
- **Numbered messages**: Sequence numbers on arrows (same as sequence diagram)
- **Spatial layout**: Related objects positioned near each other
- **Bidirectional arrows**: Show mutual communication

### What It Shows Well
✅ **Relationships**: See all connections between objects  
✅ **Communication patterns**: Identify heavily-connected objects  
✅ **Dependencies**: Spot tight coupling  
✅ **Data flow**: Follow information between components  
✅ **Message aggregation**: Multiple messages on one link  

### What It Doesn't Show Well
❌ **Precise timing**: Sequence numbers indicate order but not visual flow  
❌ **Nested calls**: Call stack depth not obvious  
❌ **Control flow**: Loops and conditionals less clear  

### Example: Refactoring Review

**Scenario**: Want to swap out Gurobi for a different LP solver

**How to use the communication diagram:**
1. Find `GrbLPSolverImp` object
2. See it has links to: `Scenario`, `Reaction`, `Metabolite`, `LPSolver`
3. Check messages: reads data, formulates problem, returns results
4. Compare with `HighsLPSolverImp` and `LpsLPSolverImp` objects
5. Verify they have identical communication patterns
6. Confirm the interface (`LPSolverImp`) is the only dependency

You can see **what needs to stay compatible** when changing implementations.

---

## 3. Communication Diagram Compact (reftest_communication_compact.puml)

### Purpose
Shows **system architecture** and **high-level organization**.

### Visual Characteristics
- **Packages**: Components grouped into logical layers
- **Inheritance**: `<|--` notation shows interface implementation
- **Simplified messages**: Aggregated descriptions instead of all details
- **Layer hierarchy**: Top-level system down to implementations
- **Design patterns**: Annotations identify architectural patterns

### What It Shows Well
✅ **Architecture layers**: See system decomposition  
✅ **Design patterns**: Strategy, Facade, Composite patterns visible  
✅ **Module boundaries**: Package groupings show separation of concerns  
✅ **Inheritance hierarchies**: Interface/implementation relationships  
✅ **Key dependencies**: Major communication paths emphasized  

### What It Doesn't Show Well
❌ **Detailed interactions**: Message content simplified  
❌ **All components**: Minor utility objects omitted  
❌ **Precise sequence**: Numbers present but not emphasized  

### Example: Architecture Presentation

**Scenario**: Explaining PLUM architecture to new developers

**How to use the compact communication diagram:**
1. Start at top: **Test Layer** (User → reftest.sh)
2. Move down: **Main Executable** (plum main orchestrates)
3. Point out: **Configuration** layer (Opts, Config, Params)
4. Explain: **Data Model** layer (Scenario, Reaction, Metabolite)
5. Highlight: **Solver Framework** (Strategy pattern for algorithms)
6. Show: **LP/MILP Implementation** (Strategy pattern for backends)
7. End with: **Results** layer (Solution output)

You can explain the **entire system** in 5 minutes with one diagram.

---

## Detailed Comparison

### Message Representation

**Sequence Diagram:**
```
PlumMain -> Scenario : read_data(data_fn)
Scenario -> Metabolite : create metabolites
Metabolite --> Scenario : metabolite objects
Scenario --> PlumMain : data loaded
```
**Shows**: Exact call sequence, activation bars, return values

**Communication Diagram:**
```
PlumMain --> Scenario : "22: read_data(data_fn)
                         23: read_flux(flux_fn)
                         24: read_react_cost(cost_fn)
                         25: read_supply_demand(sd_fn)
                         26: finalise(params)"
```
**Shows**: All related messages on one link with sequence numbers

**Communication Compact:**
```
PlumMain -down-> Scenario : "22-26: load data"
```
**Shows**: Aggregated purpose of message group

---

### Object Lifetime

**Sequence Diagram:**
- Shows when objects are `activate`d and `deactivate`d
- Vertical lifelines indicate existence throughout scenario
- Activation bars show processing periods
- Clear creation and destruction points

**Communication Diagrams:**
- Objects shown as static entities
- Lifetime not explicitly represented
- Focus on "what exists and communicates"
- Creation indicated by message (e.g., "create solution")

---

### Control Flow

**Sequence Diagram:**
```plantuml
alt Solver type selection
    PlumMain -> ConcreteSolver: new LPSolver(...)
else INT
    PlumMain -> ConcreteSolver: new IntSolver(...)
else INCR
    PlumMain -> ConcreteSolver: new IncrSolver(...)
end

loop For each experiment
    ConcreteSolver -> SolverImp: solve()
end
```
**Shows**: Explicit control structures with visual blocks

**Communication Diagrams:**
```
PlumMain --> ConcreteSolver : "45: create solver
                               (LPSolver | IntSolver | IncrSolver)
                               46: solve()"
```
**Shows**: Alternative paths described in message text

---

## When to Use Each Diagram

### Use Sequence Diagram When:

✅ **Debugging a specific issue**
   - Need to trace exact execution path
   - Want to see where data transforms
   - Looking for missing or extra calls

✅ **Understanding algorithm flow**
   - Need to see temporal dependencies
   - Want to understand control flow
   - Studying nested call patterns

✅ **Documenting a process**
   - Explaining how a feature works step-by-step
   - Writing technical specifications
   - Creating detailed design docs

✅ **Code review for logic**
   - Verifying correct call sequence
   - Checking error handling paths
   - Ensuring proper initialization order

**Example questions answered:**
- "What happens when plum is invoked with -v INT?"
- "Where does the biomass flux get calculated?"
- "In what order are files read?"
- "When does the LP solver get instantiated?"

---

### Use Communication Diagram When:

✅ **Refactoring planning**
   - Identifying coupled components
   - Finding candidates for extraction
   - Understanding impact of changes

✅ **Dependency analysis**
   - Seeing what depends on what
   - Finding circular dependencies
   - Planning module boundaries

✅ **Design review**
   - Evaluating architectural decisions
   - Checking interface contracts
   - Verifying separation of concerns

✅ **Integration work**
   - Understanding component interfaces
   - Seeing communication patterns
   - Planning API changes

**Example questions answered:**
- "What objects does Scenario communicate with?"
- "If I change Reaction, what's affected?"
- "How many components depend on Params?"
- "What's the interface between PlumMain and solvers?"

---

### Use Communication Compact Diagram When:

✅ **Presentations and talks**
   - Introducing system to audience
   - High-level architecture overview
   - Time-constrained explanations

✅ **Onboarding new developers**
   - First introduction to codebase
   - Understanding major components
   - Seeing layer organization

✅ **Documentation and reports**
   - Architecture section of papers
   - Technical reports
   - Project documentation

✅ **Strategic planning**
   - Identifying major subsystems
   - Planning architectural changes
   - Evaluating technology choices

**Example questions answered:**
- "What are the major components of PLUM?"
- "How is the system organized?"
- "What design patterns are used?"
- "Where would feature X fit in the architecture?"

---

## Practical Workflow

### Starting a New Task

1. **First**: Look at **Communication Compact**
   - Get oriented: Where does this feature fit?
   - Identify relevant layers/packages
   - Understand high-level flow

2. **Second**: Study **Communication Diagram**
   - See detailed component interactions
   - Identify which objects you'll modify
   - Understand dependencies

3. **Third**: Trace **Sequence Diagram**
   - Follow exact execution path
   - Understand timing and order
   - See where to add your code

### Debugging a Problem

1. **First**: Use **Sequence Diagram**
   - Trace the problematic execution path
   - Find where behavior deviates
   - Identify the failing call

2. **Second**: Check **Communication Diagram**
   - See what else might affect that component
   - Look for unexpected dependencies
   - Identify side effects

3. **Third**: Review **Communication Compact**
   - Understand broader context
   - See if problem is architectural
   - Plan fix at appropriate layer

### Refactoring

1. **First**: Analyze **Communication Diagram**
   - Identify tight coupling
   - Find refactoring candidates
   - Plan interface changes

2. **Second**: Check **Communication Compact**
   - Ensure change fits architecture
   - Verify layer boundaries maintained
   - Confirm pattern consistency

3. **Third**: Verify **Sequence Diagram**
   - Ensure execution flow preserved
   - Check all paths updated
   - Validate temporal contracts

---

## Diagram Maintenance

### When to Update

Update **all three** diagrams when:
- ✅ Adding new major component
- ✅ Changing component relationships
- ✅ Modifying core execution flow
- ✅ Refactoring architecture

Update **sequence and communication** when:
- ✅ Adding new method calls
- ✅ Changing data flow
- ✅ Modifying control flow

Update **only sequence** when:
- ✅ Changing method parameters
- ✅ Adjusting timing/ordering
- ✅ Adding error handling

### Consistency Rules

1. **Message numbers must match**
   - Sequence diagram defines canonical numbering
   - Communication diagrams reference same numbers
   - Compact diagram groups numbers logically

2. **Object names must be consistent**
   - Use same names across all diagrams
   - Match actual class names in code
   - Use colons for instances (`:Scenario`)

3. **Abstraction levels must align**
   - Sequence: Detailed method calls
   - Communication: Method groups
   - Compact: Logical operations

---

## Advanced Tips

### Reading Sequence Diagrams

**Tip 1**: Follow activation bars vertically
- Shows call stack depth
- Nested bars = nested calls
- Long bars = complex operations

**Tip 2**: Look for return arrows
- Dashed arrows = return values
- Short returns = simple data
- Long path returns = propagated results

**Tip 3**: Study alt/loop blocks
- Shows conditional paths
- Identifies variation points
- Reveals algorithmic structure

### Reading Communication Diagrams

**Tip 1**: Count arrows per object
- Many arrows = hub/facade
- Few arrows = leaf/utility
- Bidirectional = peer communication

**Tip 2**: Follow message numbers
- Low numbers = initialization
- Middle numbers = main work
- High numbers = cleanup/output

**Tip 3**: Identify star patterns
- Center object = orchestrator
- Spokes = delegated responsibilities
- Pattern = Facade or Mediator

### Reading Compact Communication Diagrams

**Tip 1**: Observe package organization
- Top packages = external interface
- Middle packages = business logic
- Bottom packages = infrastructure

**Tip 2**: Look for inheritance (triangles)
- Multiple implementations = Strategy
- Interface + implementations = polymorphism
- Shows extension points

**Tip 3**: Trace layer transitions
- Clear boundaries = good separation
- Cross-cutting links = coupling
- Skipped layers = architectural smell

---

## Tools and Rendering

### PlantUML Commands

```bash
# Generate all diagrams
cd doc/uml
plantuml *.puml

# Generate specific formats
plantuml -tsvg *.puml    # SVG for web
plantuml -tpng *.puml    # PNG for docs
plantuml -tpdf *.puml    # PDF for reports

# Generate with custom DPI
plantuml -DPLANTUML_LIMIT_SIZE=8192 *.puml
```

### Viewing Recommendations

**Sequence Diagram:**
- Large vertical monitor or rotated display
- Zoom to 100-150% for readability
- Use PDF for high-quality printing

**Communication Diagram:**
- Standard landscape orientation
- SVG format for pan/zoom
- Good for interactive exploration

**Communication Compact:**
- Fits well on single slide
- PNG format for presentations
- Print-friendly at letter size

---

## Conclusion

All three diagrams **complement each other**:

- **Sequence**: "What happens and when?"
- **Communication**: "Who talks to whom?"
- **Compact**: "How is it organized?"

Use them together for **complete understanding** of PLUM's architecture and execution.

**Quick decision tree:**
```
Are you debugging or tracing execution?
├─ YES → Use Sequence Diagram
└─ NO → Are you reviewing architecture?
    ├─ YES → Start with Compact, then Communication
    └─ NO → Are you planning changes?
        ├─ YES → Use Communication Diagram
        └─ NO → Want overview? → Use Compact
```

---

## Additional Resources

- **This README**: `doc/uml/README.md` (generation instructions)
- **Creation story**: `doc/uml/SUMMARY.md` (feasibility analysis)
- **PlantUML docs**: https://plantuml.com/
- **UML spec**: https://www.omg.org/spec/UML/
- **PLUM architecture**: `CLAUDE.md`

---

**Created**: 2026-06-15  
**Author**: Claude Code (Anthropic)  
**Project**: PLUM Metabolic Gap-Filling Solver
