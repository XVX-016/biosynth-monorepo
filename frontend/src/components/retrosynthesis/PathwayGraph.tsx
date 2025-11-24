import React from 'react'
import { useLabStore } from '../../../store/labStore'

interface PathwayStep {
  molecule: any
  step: number
  reaction?: any
  precursors?: any[]
  is_starting?: boolean
}

interface PathwayGraphProps {
  pathway: {
    steps: PathwayStep[]
    score?: number
    score_breakdown?: any
  }
  onStepClick?: (step: PathwayStep) => void
}

export default function PathwayGraph({ pathway, onStepClick }: PathwayGraphProps) {
  const steps = pathway.steps || []

  return (
    <div className="w-full space-y-3">
      <div className="text-xs font-semibold text-gray-700 mb-2">
        Pathway Score: {pathway.score?.toFixed(2) || 'N/A'}
      </div>
      
      {pathway.score_breakdown && (
        <div className="text-xs text-gray-600 mb-3 space-y-1">
          <div>Step Score: {pathway.score_breakdown.step_score?.toFixed(2)}</div>
          <div>Reagent Score: {pathway.score_breakdown.reagent_score?.toFixed(2)}</div>
          <div>Energy Score: {pathway.score_breakdown.energy_score?.toFixed(2)}</div>
          <div>Yield Score: {pathway.score_breakdown.yield_score?.toFixed(2)}</div>
        </div>
      )}

      <div className="space-y-2">
        {steps.map((step, index) => (
          <div key={index} className="border border-gray-200 rounded p-2">
            <div className="flex items-center justify-between mb-1">
              <span className="text-xs font-semibold text-gray-700">
                Step {step.step}
                {step.is_starting && <span className="text-green-600 ml-1">(Starting Material)</span>}
              </span>
              {step.molecule && (
                <span className="text-xs text-gray-500">
                  {step.molecule.atoms?.length || 0} atoms
                </span>
              )}
            </div>
            
            {step.reaction && (
              <div className="text-xs text-gray-600 mb-1">
                Reaction: {step.reaction.name || step.reaction.id}
              </div>
            )}
            
            {step.precursors && step.precursors.length > 0 && (
              <div className="text-xs text-gray-600 mb-1">
                Precursors: {step.precursors.length}
              </div>
            )}
            
            <button
              onClick={() => {
                if (step.molecule && onStepClick) {
                  onStepClick(step)
                }
              }}
              className="mt-1 px-2 py-1 text-xs bg-blue-600 text-white rounded hover:bg-blue-700"
            >
              Load in Viewer
            </button>
          </div>
        ))}
      </div>
    </div>
  )
}

