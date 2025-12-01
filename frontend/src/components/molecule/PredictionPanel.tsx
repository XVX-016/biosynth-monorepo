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
import type { PredictionResult, AttentionData } from '@/lib/molecule/prediction'

interface PredictionPanelProps {
  molecule: Molecule
  onAttentionUpdate?: (attention: AttentionData | null) => void
  className?: string
}

export function PredictionPanel({
  molecule,
  onAttentionUpdate,
  className = '',
}: PredictionPanelProps) {
  const [predictions, setPredictions] = useState<PredictionResult | null>(null)
  const [loading, setLoading] = useState(false)
  const [error, setError] = useState<string | null>(null)
  const [selectedProperties, setSelectedProperties] = useState<string[]>(['logP', 'solubility', 'toxicity'])
  const abortControllerRef = useRef<AbortController | null>(null)

  // Predict when molecule changes
  useEffect(() => {
    if (molecule.isEmpty()) {
      setPredictions(null)
      setError(null)
      onAttentionUpdate?.(null)
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
          undefined, // Use default model
          selectedProperties
        )

        if (!abortControllerRef.current?.signal.aborted) {
          setPredictions(result)
          onAttentionUpdate?.(result.attention || null)
        }
      } catch (err) {
        if (!abortControllerRef.current?.signal.aborted) {
          const errorMessage = err instanceof Error ? err.message : 'Prediction failed'
          setError(errorMessage)
          setPredictions(null)
          onAttentionUpdate?.(null)
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
  }, [molecule, selectedProperties, onAttentionUpdate])

  // Cleanup on unmount
  useEffect(() => {
    return () => {
      predictionService.cancel()
      if (abortControllerRef.current) {
        abortControllerRef.current.abort()
      }
    }
  }, [])

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

          {predictions.attention && (
            <div className="mt-3 pt-3 border-t border-gray-200">
              <div className="text-xs text-gray-500">
                Attention map available ({predictions.attention.edgeAttentions.length} edges)
              </div>
            </div>
          )}
        </div>
      )}

      {!predictions && !loading && !error && (
        <div className="text-xs text-gray-500">Waiting for molecule...</div>
      )}
    </div>
  )
}

