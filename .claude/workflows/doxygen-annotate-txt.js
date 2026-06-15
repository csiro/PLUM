export const meta = {
  name: 'doxygen-annotate',
  description: 'Add Doxygen documentation to all C/C++ source and header files',
  phases: [
    { title: 'Analyze', detail: 'Read and understand each file structure' },
    { title: 'Document', detail: 'Add Doxygen comments to files' },
  ],
}

const files = [
  // Header files
  './include/mosh/builddate.h',
  './include/mosh/constants.h',
  './include/mosh/criticalpath.h',
  './include/mosh/critpathnode.h',
  './include/mosh/dualvals.h',
  './include/mosh/dummysolver.h',
  './include/mosh/experiment.h',
  './include/mosh/gapsolver.h',
  './include/mosh/grblpsolverimp.h',
  './include/mosh/highslpsolverimp.h',
  './include/mosh/highsutil.h',
  './include/mosh/incrsolver.h',
  './include/mosh/intsolver.h',
  './include/mosh/lnsmxsolver.h',
  './include/mosh/lpslpsolverimp.h',
  './include/mosh/lpsolver.h',
  './include/mosh/lpsolverimp.h',
  './include/mosh/lpsutil.h',
  './include/mosh/mathheur.h',
  './include/mosh/metabolite.h',
  './include/mosh/multisol.h',
  './include/mosh/params.h',
  './include/mosh/pathfinder.h',
  './include/mosh/path.h',
  './include/mosh/pathsol.h',
  './include/mosh/reaction.h',
  './include/mosh/scenario.h',
  './include/mosh/solution.h',
  './include/mosh/supplyresid.h',
  // Source files
  './src/criticalpath.cpp',
  './src/critpathnode.cpp',
  './src/dummysolver.cpp',
  './src/grblpsolverimp.cpp',
  './src/highslpsolverimp.cpp',
  './src/highsutil.cpp',
  './src/incrsolver.cpp',
  './src/intsolver.cpp',
  './src/lnsmxsolver.cpp',
  './src/lpslpsolverimp.cpp',
  './src/lpsolver.cpp',
  './src/lpsutil.cpp',
  './src/mathheur.cpp',
  './src/merge.cpp',
  './src/mergeplum.cpp',
  './src/multisol.cpp',
  './src/overlap.cpp',
  './src/params.cpp',
  './src/pathfinder.cpp',
  './src/plumchk.cpp',
  './src/plumcmp.cpp',
  './src/plum.cpp',
  './src/plummerge.cpp',
  './src/plummx.cpp',
  './src/plumseq.cpp',
  './src/plumsp.cpp',
  './src/scenario.cpp',
  './src/solution.cpp',
  './src/strex.cpp',
]

const DOXYGEN_SCHEMA = {
  type: 'object',
  required: ['annotations'],
  properties: {
    annotations: {
      type: 'array',
      description: 'List of Doxygen annotations to add',
      items: {
        type: 'object',
        required: ['line_number', 'doxygen_comment', 'target_type'],
        properties: {
          line_number: {
            type: 'number',
            description: 'Line number where the Doxygen comment should be inserted (before the declaration)'
          },
          doxygen_comment: {
            type: 'string',
            description: 'Complete Doxygen comment block to insert'
          },
          target_type: {
            type: 'string',
            description: 'Type of element being documented: file/class/function/variable/enum/typedef'
          },
          target_name: {
            type: 'string',
            description: 'Name of the element being documented'
          }
        }
      }
    }
  }
}

phase('Analyze')
log(`Analyzing ${files.length} files to add Doxygen documentation`)

const analyses = await pipeline(
  files,
  (file) => agent(`Read the file ${file} and analyze its structure.

For this file, identify ALL elements that need Doxygen documentation:
- File header (should be at line 1)
- Classes and structs
- All public, protected, and private member functions
- Member variables
- Free functions
- Enums and typedefs
- Namespaces

For each element, provide:
1. The line number where it appears
2. A complete Doxygen comment block using appropriate format:
   - File: /** @file */ with @brief description
   - Class/struct: /** @class/@struct */ with @brief and detailed description
   - Function: /** @brief */ with @param and @return tags
   - Variable: /**< Brief description */ or /** @brief */
   - Enum: /** @enum */ with @brief for enum and values
   
Consider the context of metabolic gap-filling and flux balance analysis. Use domain terminology appropriately.

Return a structured list of annotations sorted by line number (ascending).`, {
    label: file.split('/').pop(),
    phase: 'Analyze',
    schema: DOXYGEN_SCHEMA
  })
)

phase('Document')
log(`Adding Doxygen comments to ${analyses.filter(Boolean).length} files`)

const results = await pipeline(
  analyses.filter(Boolean).map((analysis, idx) => ({ analysis, file: files[idx] })),
  ({ analysis, file }) => agent(`Add Doxygen documentation to the file ${file}.

You have the following annotations to add:
${JSON.stringify(analysis.annotations, null, 2)}

Your task:
1. Read the current file
2. Add each Doxygen comment at the specified line number
3. Process annotations from highest line number to lowest (so line numbers don't shift)
4. Ensure proper formatting and indentation
5. Make sure the Doxygen comments follow standard conventions

For each annotation, insert the doxygen_comment text at the specified line_number.

After adding all annotations, return a brief summary of what was documented.`, {
    label: `doc:${file.split('/').pop()}`,
    phase: 'Document'
  })
)

const successful = results.filter(Boolean)
log(`Successfully documented ${successful.length} files`)

return {
  total_files: files.length,
  documented: successful.length,
  summary: `Added Doxygen documentation to ${successful.length} C/C++ files`
}
