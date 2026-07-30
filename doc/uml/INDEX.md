# PLUM UML Diagrams - Complete Index

This directory contains comprehensive UML diagrams documenting PLUM's reference test execution flow from bash script to C++ class level.

---

## 📊 Diagrams Available

### Behavioral Diagrams (Execution Flow)

#### 1. Sequence Diagram
**File**: `reftest_sequence.puml` (315 lines)  
**Generates**: `reftest_sequence.png`, `reftest_sequence.svg`

**Focus**: Temporal execution flow  
**Shows**: Chronological order of method calls, control flow, activation periods  
**Use for**: Debugging, tracing execution, understanding algorithms  

#### 2. Communication Diagram (Detailed)
**File**: `reftest_communication.puml` (239 lines)  
**Generates**: `reftest_communication.png`, `reftest_communication.svg`

**Focus**: Object relationships and interactions  
**Shows**: Communication links, message flow, dependencies  
**Use for**: Design review, refactoring, dependency analysis  

#### 3. Communication Diagram (Compact)
**File**: `reftest_communication_compact.puml` (188 lines)  
**Generates**: `reftest_communication_compact.png`, `reftest_communication_compact.svg`

**Focus**: High-level architecture  
**Shows**: Package organization, design patterns, layer hierarchy  
**Use for**: Presentations, onboarding, architecture overview  

---

### Structural Diagrams (Class Architecture)

#### 4. Class Diagram: Solver Hierarchy
**File**: `class_solver_hierarchy.puml` (343 lines)  
**Generates**: `class_solver_hierarchy.png`, `class_solver_hierarchy.svg`

**Focus**: Solver framework and LP/MILP implementations  
**Shows**: GapSolver hierarchy (LPSolver, IntSolver, IncrSolver), LPSolverImp hierarchy (Gurobi, HiGHS, lp_solve), design patterns  
**Use for**: Implementing solvers, understanding algorithms, extending backends  

#### 5. Class Diagram: Data Model
**File**: `class_data_model.puml` (385 lines)  
**Generates**: `class_data_model.png`, `class_data_model.svg`

**Focus**: Metabolic network data structures  
**Shows**: Scenario, Metabolite, Reaction, Experiment, Solution, Params, complete attributes and relationships  
**Use for**: Data structure understanding, file I/O, feature implementation  

#### 6. Class Diagram: System Overview
**File**: `class_system_overview.puml` (357 lines)  
**Generates**: `class_system_overview.png`, `class_system_overview.svg`

**Focus**: Complete system architecture with packages  
**Shows**: All major components, 6 logical packages, system data flow  
**Use for**: System overview, presentations, architectural planning  

#### 7. Class Diagram: Single-Path Solver Slice
**File**: `GapSolver-To-HighsLPSolverImp.puml`  
**Generates**: `GapSolver-To-HighsLPSolverImp.svg`

**Focus**: One concrete path — `plum -v CTS -V HIGHS`  
**Shows**: Only GapSolver → LPSolver → LPSolverImp → HighsLPSolverImp (plus Scenario, Params); all sibling solvers and backends omitted  
**Use for**: Learning or debugging a single code path without full-hierarchy noise  

---

## 📚 Documentation Files

### README.md (421 lines)
**Complete user guide for the diagrams**

Contents:
- Diagram descriptions and comparison table
- PlantUML installation instructions
- Image generation commands (PNG, SVG, PDF, etc.)
- How to read sequence and communication diagrams
- Online rendering options
- Editor integration guide
- Maintenance and contribution guidelines
- PlantUML syntax reference

**Start here if**: You want to generate the diagrams or understand how to use them

---

### DIAGRAM_GUIDE.md (526 lines)
**In-depth guide comparing behavioral diagram types**

Contents:
- Quick reference table
- Detailed explanation of sequence vs communication diagrams
- What each diagram shows well (and doesn't)
- Practical examples for debugging, refactoring, presenting
- When to use which diagram (decision tree)
- Message representation comparison
- Practical workflow recommendations
- Maintenance guidelines
- Advanced reading tips
- Tools and rendering recommendations

**Start here if**: You want to understand the differences between behavioral diagrams or learn how to read them effectively

---

### CLASS_DIAGRAM_GUIDE.md (567 lines)
**In-depth guide for class diagrams**

Contents:
- Quick reference table for class diagrams
- Detailed explanation of each class diagram
- UML notation guide (access modifiers, relationships, multiplicity)
- Reading strategies (top-down, bottom-up)
- Practical workflows (implementing features, debugging, refactoring)
- Comparison with behavioral diagrams
- Advanced tips (finding information, understanding design decisions)
- Extension point identification
- Tools and rendering recommendations
- Maintenance guidelines

**Start here if**: You want to understand the class diagrams in depth or learn how to use them for development

---

### SUMMARY.md (462 lines)
**Creation story and feasibility analysis**

Contents:
- Feasibility assessment: "Is this possible?" → **YES!**
- Technical approach explanation
- Source code analysis details
- Challenges and solutions
- Diagram characteristics (participants, sections, interactions)
- Value proposition (for understanding, development, documentation)
- Validation checklist
- Limitations and considerations
- Future enhancements (class diagrams, activity diagrams, etc.)
- Files created summary

**Start here if**: You want to understand how these diagrams were created or are considering creating similar diagrams for other projects

---

### INDEX.md (this file)
**Quick navigation and overview**

---

## 🚀 Quick Start

### Generate All Diagrams (no installation required)

Use the `render.sh` helper — it needs no `apt-get` and no Graphviz. It uses the
bundled PlantUML jar plus PlantUML's built-in pure-Java Smetana layout engine
(enabled by the `!pragma layout smetana` line in each `.puml`).

```bash
# Render every *.puml to SVG (output in doc/uml/svg/)
/home/kubeflow/PLUM/doc/uml/render.sh

# Render a single diagram
/home/kubeflow/PLUM/doc/uml/render.sh GapSolver-To-HighsLPSolverImp.puml

# View results
ls -lh /home/kubeflow/PLUM/doc/uml/svg/*.svg
```

If you have permission to install packages, a system-wide PlantUML also works
(`sudo apt-get install plantuml graphviz`; then `plantuml -tsvg *.puml`). The
Smetana pragma is still honored and makes Graphviz optional.

### View Online (No Installation)

1. Go to: https://www.plantuml.com/plantuml/uml/
2. Copy contents of any `.puml` file
3. Paste and render

---

## 📖 Reading Guide

### First Time Viewing?

**Recommended order:**

1. **Start with System Overview Class Diagram**
   - Get 10,000-foot view of entire architecture
   - Understand major packages and components
   - See design patterns used
   - **File**: `class_system_overview.puml`

2. **Then Compact Communication Diagram**
   - See execution flow at high level
   - Understand layer interactions
   - Follow data flow
   - **File**: `reftest_communication_compact.puml`

3. **Dive into specific class diagrams**
   - **Solver Hierarchy**: If working with algorithms
   - **Data Model**: If working with data structures
   - **Files**: `class_solver_hierarchy.puml`, `class_data_model.puml`

4. **For detailed execution understanding**
   - **Communication Diagram**: Object relationships during execution
   - **Sequence Diagram**: Chronological method calls
   - **Files**: `reftest_communication.puml`, `reftest_sequence.puml`

### Working on Specific Task?

**Understanding the architecture?**
→ **System overview class diagram** (`class_system_overview.puml`)

**Implementing a new solver?**
→ **Solver hierarchy class diagram** (`class_solver_hierarchy.puml`)

**Working with data structures?**
→ **Data model class diagram** (`class_data_model.puml`)

**Debugging execution flow?**
→ **Sequence diagram** (`reftest_sequence.puml`)

**Planning a refactoring?**
→ **Class diagrams** + **Communication diagram** (`class_*.puml`, `reftest_communication.puml`)

**Giving a presentation?**
→ **System overview** + **Compact communication** (`class_system_overview.puml`, `reftest_communication_compact.puml`)

**New to PLUM?**
→ **All diagrams**, in recommended order above, plus read `CLASS_DIAGRAM_GUIDE.md` and `DIAGRAM_GUIDE.md`

---

## 📊 Diagram Statistics

### Coverage

**Components**: 20+ objects across 7 logical layers
- Test Layer: User, reftest.sh
- Main Executable: plum main()
- Configuration: Opts, Config, Params, Rand
- Data Model: Scenario, Metabolite, Reaction, Experiment
- Solver Framework: GapSolver, LPSolver, IntSolver, IncrSolver
- LP/MILP Implementation: LPSolverImp, GrbLPSolverImp, HighsLPSolverImp, LpsLPSolverImp
- Results: Solution, MultiSol

**Messages**: 86 numbered interactions
1-7: Test initialization
8-21: Main setup & configuration
22-38: Scenario loading
39-44: Solver implementation setup
45-64: Problem solving
65-81: Solution extraction
82-86: Test validation

**Execution Phases**: 11 major sections
1. Test Initialization
2. Test Execution Loop
3. PLUM Executable Invocation
4. Random Seed Setup
5. Scenario Loading
6. Solver Setup
7. Problem Solving
8. Results Processing
9. Output Generation
10. Test Validation
11. Test Completion

### File Sizes

| File | Lines | Size | Purpose |
|------|-------|------|---------|
| **Behavioral Diagrams** | | | |
| `reftest_sequence.puml` | 315 | 12KB | Sequence diagram |
| `reftest_communication.puml` | 239 | 6.9KB | Communication diagram |
| `reftest_communication_compact.puml` | 188 | 5.5KB | Compact communication |
| **Structural Diagrams** | | | |
| `class_solver_hierarchy.puml` | 343 | 11KB | Solver class diagram |
| `class_data_model.puml` | 385 | 12KB | Data model class diagram |
| `class_system_overview.puml` | 357 | 11KB | System overview class diagram |
| **Documentation** | | | |
| `README.md` | 558 | 18KB | User guide |
| `DIAGRAM_GUIDE.md` | 526 | 16KB | Behavioral diagram guide |
| `CLASS_DIAGRAM_GUIDE.md` | 567 | 18KB | Class diagram guide |
| `SUMMARY.md` | 462 | 16KB | Creation story |
| `INDEX.md` | (this file) | 14KB | Navigation |
| **Total** | **~4,800** | **~140KB** | Complete documentation |

---

## 🎯 Use Cases by Role

### Software Developer
**Daily work**: Debugging, feature development, code review
**Primary diagram**: Sequence → See exact execution flow
**Also useful**: Communication → Understand dependencies

**Workflow**:
1. Read feature request
2. Check Compact Communication → Find relevant layer
3. Study Communication → Identify affected components
4. Trace Sequence → See where to add code

---

### Software Architect
**Daily work**: Design decisions, architecture reviews, technical debt
**Primary diagram**: Communication Compact → See structure
**Also useful**: Communication → Analyze coupling

**Workflow**:
1. Identify architectural concern
2. Check Compact Communication → See current organization
3. Study Communication → Analyze dependencies
4. Plan changes maintaining patterns

---

### Technical Lead
**Daily work**: Code review, mentoring, documentation
**Primary diagram**: All three → Complete understanding
**Use for**: Teaching, reviewing, planning

**Workflow**:
1. Review PR → Trace in Sequence diagram
2. Mentor developer → Explain with Compact
3. Plan sprint → Check Communication for scope

---

### Researcher / Student
**Daily work**: Understanding metabolic solvers, FBA, gap-filling
**Primary diagram**: Compact Communication → Domain understanding
**Also useful**: Sequence → Algorithm details

**Workflow**:
1. Read PLUM paper
2. View Compact Communication → See system design
3. Study Sequence → Understand algorithm flow
4. Read code with diagrams as reference

---

### Project Manager / Stakeholder
**Daily work**: Progress tracking, risk assessment, planning
**Primary diagram**: Compact Communication → High-level view
**Use for**: Understanding scope, identifying risks

**Workflow**:
1. View Compact Communication → Understand system
2. Identify components affected by requirements
3. Estimate effort based on complexity
4. Track progress against architecture

---

## 🔄 Relationships Between Files

```
INDEX.md (you are here)
├─ Quick navigation to all resources
├─ Usage recommendations by role
└─ Points to all other files

README.md
├─ Installation and generation instructions
├─ Diagram descriptions
└─ Technical how-to guide

DIAGRAM_GUIDE.md
├─ Deep comparison of diagram types
├─ When to use which
├─ Practical workflows
└─ Advanced tips

SUMMARY.md
├─ Creation story
├─ Feasibility analysis
├─ Technical approach
└─ Future enhancements

reftest_sequence.puml
├─ Source: Sequence diagram
├─ Generates: PNG, SVG, PDF
└─ Referenced by: All documentation

reftest_communication.puml
├─ Source: Communication diagram (detailed)
├─ Generates: PNG, SVG, PDF
└─ Uses same message numbers as sequence

reftest_communication_compact.puml
├─ Source: Communication diagram (compact)
├─ Generates: PNG, SVG, PDF
└─ Simplified version of communication
```

---

## 🛠️ Maintenance

### Keeping Diagrams Updated

**When code changes, update diagrams:**

| Change Type | Update Required |
|-------------|----------------|
| New class added | All three diagrams |
| Class removed | All three diagrams |
| Method signature changed | Sequence only |
| New relationship | Communication diagrams |
| Architecture change | All three diagrams |
| Algorithm refinement | Sequence + Communication |
| Bug fix (no structure change) | No update needed |

**Validation checklist:**
- [ ] Message numbers consistent across diagrams
- [ ] Object names match code (`:ClassName` format)
- [ ] All major components represented
- [ ] Inheritance relationships correct
- [ ] PlantUML syntax validates (test with `plantuml -checkonly`)
- [ ] Generated images render correctly

---

## 📦 Integration with PLUM Documentation

These diagrams are part of PLUM's comprehensive documentation:

```
PLUM/
├── CLAUDE.md                          # Project overview and build guide
├── doc/
│   ├── overview/overview.md          # Project background and methodology
│   ├── fileformat.md                 # PLD file format specification
│   └── uml/                          # ⭐ You are here
│       ├── INDEX.md                  # This file
│       ├── README.md                 # User guide
│       ├── DIAGRAM_GUIDE.md          # Comparison guide
│       ├── SUMMARY.md                # Creation story
│       ├── reftest_sequence.puml     # Sequence diagram
│       ├── reftest_communication.puml # Communication diagram
│       └── reftest_communication_compact.puml # Compact communication
├── include/mosh/                     # Header files (Doxygen documented)
├── src/                              # Source files (Doxygen documented)
├── util/                             # Python utilities (numpy-style docstrings)
└── reftest/
    └── reftest.sh                    # Test harness (bash documentation)
```

**See also**:
- `CLAUDE.md` → Build instructions, architecture overview
- `doc/overview/overview.md` → Project background, methodology
- Header files → Doxygen API documentation
- Python scripts → Numpy-style docstrings

---

## 🌐 External Resources

### PlantUML
- **Website**: https://plantuml.com/
- **Sequence diagrams**: https://plantuml.com/sequence-diagram
- **Communication diagrams**: PlantUML uses object diagram syntax
- **Online editor**: https://www.plantuml.com/plantuml/uml/

### UML Standard
- **OMG UML Spec**: https://www.omg.org/spec/UML/
- **Wikipedia**: https://en.wikipedia.org/wiki/Unified_Modeling_Language
- **Tutorials**: https://www.uml-diagrams.org/

### PLUM Project
- **GitHub**: (add repository URL if public)
- **CSIRO**: https://www.csiro.au/
- **FBA Background**: (add relevant publications)

---

## 💬 Questions or Issues?

### For diagram generation issues:
→ Check `README.md` section "Generating Diagrams"

### For understanding what diagrams show:
→ Read `DIAGRAM_GUIDE.md` section "When to Use Each Diagram"

### For technical details on creation:
→ See `SUMMARY.md` section "Technical Approach"

### For PlantUML syntax:
→ Visit https://plantuml.com/

### For PLUM architecture questions:
→ See `/home/kubeflow/PLUM/CLAUDE.md`

---

## 📝 Version History

**v1.0 (2026-06-15)**
- Initial creation of all three diagram types
- Complete documentation suite
- Integration with PLUM project documentation
- Comprehensive user guides and comparison analysis

**Future additions** (see SUMMARY.md):
- Class diagrams (solver hierarchy, data model)
- Activity diagrams (reachability analysis, incremental solver)
- Component diagrams (system architecture)
- State diagrams (solver states)
- Deployment diagrams (dependencies and environments)

---

## 📄 License

These diagrams and documentation are part of the PLUM project.
(Add appropriate license information)

---

## ✨ Credits

**Created**: 2026-06-15  
**Author**: Claude Code (Anthropic)  
**Project**: PLUM - Metabolic Gap-Filling Solver  
**Organization**: CSIRO  

**Acknowledgments**:
- PLUM development team for well-documented codebase
- Doxygen documentation in all C++ files
- Numpy-style docstrings in Python utilities
- Comprehensive CLAUDE.md for architecture context

---

**Quick Navigation**:
- [📖 User Guide](README.md)
- [📊 Diagram Comparison](DIAGRAM_GUIDE.md)
- [🔍 Creation Story](SUMMARY.md)
- [🏗️ Project Documentation](../../CLAUDE.md)

**Generate diagrams**: `doc/uml/render.sh` (no install needed; bundled jar + Smetana)  
**View online**: https://www.plantuml.com/plantuml/uml/

---

*Last updated: 2026-06-15*
