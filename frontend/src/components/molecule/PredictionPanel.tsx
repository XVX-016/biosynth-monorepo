/**
 * PredictionPanel - ML Prediction UI Component
 * 
 * Phase 7: ML Prediction Pipeline Fix
 * 
 * Features:
 * - Debounced predictions (300ms)
 * - Real-time updates as molecule changes
 * - Attention map visualization
 * - Error handling
 * - Loading states
 */

import React, { useState, useEffect, useCallback, useRef } from 'react'
import { Molecule, predictionService } from '@/lib/molecule'
import type { PredictionResult } from '@/lib/molecule/prediction'
import { useEditorContext } from './EditorContext'

interface PredictionPanelProps {
  molecule: Molecule
  className?: string
}

const AVAILABLE_MODELS = [
  { id: 'propnet', label: 'PropNet (GAT)' },
  { id: 'chem_gpt', label: 'ChemGPT (beta)' },
]

export function PredictionPanel({
  molecule,
  className = '',
}: PredictionPanelProps) {
  const [predictions, setPredictions] = useState<PredictionResult | null>(null)
  const [loading, setLoading] = useState(false)
  const [error, setError] = useState<string | null>(null)
  const [selectedProperties, setSelectedProperties] = useState<string[]>(['logP', 'solubility', 'toxicity'])
  const abortControllerRef = useRef<AbortController | null>(null)
  const {
    setAttentionData,
    attentionOverlayEnabled,
    setAttentionOverlayEnabled,
    attentionMap,
    selectedModelId,
    setSelectedModelId,
  } = useEditorContext()

  // Predict when molecule changes
  useEffect(() => {
    if (molecule.isEmpty()) {
      setPredictions(null)
      setError(null)
      setAttentionData(null)
      return
    }

    // Cancel previous request
    if (abortControllerRef.current) {
      abortControllerRef.current.abort()
    }

    // Create new abort controller
    abortControllerRef.current = new AbortController()

    // Debounced prediction
    const timeoutId = setTimeout(async () => {
      setLoading(true)
      setError(null)

      try {
        const result = await predictionService.predictWithAttention(
          molecule,
          selectedModelId,
          selectedProperties
        )

        if (!abortControllerRef.current?.signal.aborted) {
          setPredictions(result)
          setAttentionData(result.attention || null)
        }
      } catch (err) {
        if (!abortControllerRef.current?.signal.aborted) {
          const errorMessage = err instanceof Error ? err.message : 'Prediction failed'
          setError(errorMessage)
          setPredictions(null)
          setAttentionData(null)
        }
      } finally {
        if (!abortControllerRef.current?.signal.aborted) {
          setLoading(false)
        }
      }
    }, 300) // Debounce matches PredictionService

    return () => {
      clearTimeout(timeoutId)
      if (abortControllerRef.current) {
        abortControllerRef.current.abort()
      }
    }
  }, [molecule, selectedProperties, selectedModelId, setAttentionData])

  // Cleanup on unmount
  useEffect(() => {
    return () => {
      predictionService.cancel()
      if (abortControllerRef.current) {
        abortControllerRef.current.abort()
      }
      setAttentionData(null)
    }
  }, [setAttentionData])

  const formatValue = (value: number): string => {
    if (value === undefined || value === null) return 'N/A'
    return value.toFixed(3)
  }

  if (molecule.isEmpty()) {
    return (
      <div className={`bg-white rounded-lg border border-gray-200 p-4 ${className}`}>
        <h3 className="text-sm font-semibold text-gray-700 mb-2">ML Predictions</h3>
        <p className="text-xs text-gray-500">Add atoms to run predictions</p>
      </div>
    )
  }

  return (
    <div className={`bg-white rounded-lg border border-gray-200 p-4 ${className}`}>
      <div className="flex items-center justify-between mb-3">
        <h3 className="text-sm font-semibold text-gray-700">ML Predictions</h3>
        {loading && (
          <div className="flex items-center gap-2 text-xs text-gray-500">
            <div className="w-3 h-3 border-2 border-blue-500 border-t-transparent rounded-full animate-spin" />
            <span>Predicting...</span>
          </div>
        )}
      </div>

      <div className="flex flex-col gap-2 mb-3">
        <label className="text-xs text-gray-500 font-medium">Model</label>
        <select
          value={selectedModelId}
          onChange={(e) => setSelectedModelId(e.target.value)}
          className="w-full text-xs border border-gray-300 rounded px-2 py-1.5 focus:outline-none focus:ring-1 focus:ring-blue-500"
        >
          {AVAILABLE_MODELS.map((model) => (
            <option key={model.id} value={model.id}>
              {model.label}
            </option>
          ))}
        </select>
      </div>

      <div className="flex items-center justify-between mb-4">
        <span className="text-xs text-gray-600">Attention Overlay</span>
        <label className="inline-flex items-center cursor-pointer">
          <span className="mr-2 text-[11px] text-gray-500">{attentionOverlayEnabled ? 'On' : 'Off'}</span>
          <input
            type="checkbox"
            className="sr-only"
            checked={attentionOverlayEnabled}
            onChange={(e) => setAttentionOverlayEnabled(e.target.checked)}
          />
          <div
            className={`w-10 h-5 flex items-center rounded-full p-0.5 transition-colors ${
              attentionOverlayEnabled ? 'bg-blue-500/70' : 'bg-gray-300'
            }`}
          >
            <div
              className={`bg-white w-4 h-4 rounded-full shadow transform duration-200 ${
                attentionOverlayEnabled ? 'translate-x-5' : ''
              }`}
            />
          </div>
        </label>
      </div>

      {error && (
        <div className="mb-3 p-2 bg-red-50 border border-red-200 rounded text-xs text-red-700">
          {error}
        </div>
      )}

      {predictions && !loading && (
        <div className="space-y-2">
          {selectedProperties.map((prop) => {
            const value = predictions.properties[prop]
            if (value === undefined) return null

            return (
              <div key={prop} className="flex items-center justify-between text-xs">
                <span className="text-gray-600 capitalize">{prop}:</span>
                <span className="font-mono font-semibold text-gray-900">
                  {formatValue(value)}
                </span>
              </div>
            )
          })}

          {attentionMap && attentionOverlayEnabled && (
            <div className="mt-3 pt-3 border-t border-gray-200">
              <div className="text-xs text-gray-500 flex items-center justify-between mb-1">
                <span>Attention intensity</span>
                <span>
                  {attentionMap.minValue.toFixed(2)} – {attentionMap.maxValue.toFixed(2)}
                </span>
              </div>
              <div
                className="h-2 rounded-full"
                style={{
                  background:
                    'linear-gradient(90deg, rgb(0,0,255) 0%, rgb(0,255,255) 33%, rgb(255,255,0) 66%, rgb(255,0,0) 100%)',
                }}
              />
            </div>
          )}
        </div>
      )}

      {!predictions && !loading && !error && (
        <div className="text-xs text-gray-500">Waiting for molecule...</div>
      )}
      </div>
    </ErrorBoundary>
  )
}

