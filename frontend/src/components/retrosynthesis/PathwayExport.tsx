import React from 'react'

interface PathwayExportProps {
  pathway: any
  onExport?: (format: 'json' | 'csv' | 'png') => void
}

export default function PathwayExport({ pathway, onExport }: PathwayExportProps) {
  const handleExport = (format: 'json' | 'csv' | 'png') => {
    if (onExport) {
      onExport(format)
      return
    }

    // Default export implementation
    if (format === 'json') {
      const dataStr = JSON.stringify(pathway, null, 2)
      const blob = new Blob([dataStr], { type: 'application/json' })
      const url = URL.createObjectURL(blob)
      const link = document.createElement('a')
      link.href = url
      link.download = `pathway_${Date.now()}.json`
      link.click()
      URL.revokeObjectURL(url)
    } else if (format === 'csv') {
      const steps = pathway.steps || []
      const csvRows = [
        ['Step', 'Molecule', 'Reaction', 'Precursors', 'Reagents'].join(',')
      ]
      
      steps.forEach((step: any, index: number) => {
        const reaction = step.reaction?.name || step.reaction?.id || 'N/A'
        const precursors = step.precursors?.map((p: any) => p.reagent).join('; ') || 'N/A'
        csvRows.push([
          index + 1,
          `${step.molecule?.atoms?.length || 0} atoms`,
          reaction,
          step.precursors?.length || 0,
          precursors
        ].join(','))
      })
      
      const csv = csvRows.join('\n')
      const blob = new Blob([csv], { type: 'text/csv' })
      const url = URL.createObjectURL(blob)
      const link = document.createElement('a')
      link.href = url
      link.download = `pathway_${Date.now()}.csv`
      link.click()
      URL.revokeObjectURL(url)
    }
  }

  return (
    <div className="flex gap-2">
      <button
        onClick={() => handleExport('json')}
        className="px-2 py-1 text-xs bg-gray-100 hover:bg-gray-200 rounded"
      >
        Export JSON
      </button>
      <button
        onClick={() => handleExport('csv')}
        className="px-2 py-1 text-xs bg-gray-100 hover:bg-gray-200 rounded"
      >
        Export CSV
      </button>
    </div>
  )
}

