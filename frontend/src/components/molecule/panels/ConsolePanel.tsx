/**
 * ConsolePanel - Validation messages console
 * 
 * Phase 10: Lab Page UI Rebuild
 */

import React from 'react'
import { useEditorContext } from '../EditorContext'

export function ConsolePanel() {
  const { validationResult } = useEditorContext()

  if (!validationResult) {
    return (
      <div className="bg-gray-50 border border-gray-200 rounded-lg p-3">
        <h3 className="text-xs font-semibold text-gray-700 mb-2">Console</h3>
        <p className="text-xs text-gray-500">No validation messages</p>
      </div>
    )
  }

  const hasErrors = validationResult.errors.length > 0
  const hasWarnings = validationResult.warnings.length > 0

  return (
    <div className="bg-gray-50 border border-gray-200 rounded-lg p-3 max-h-32 overflow-y-auto">
      <h3 className="text-xs font-semibold text-gray-700 mb-2">Console</h3>
      
      {validationResult.valid && !hasWarnings && (
        <p className="text-xs text-green-600">✓ Molecule is valid</p>
      )}

      {hasErrors && (
        <div className="space-y-1">
          {validationResult.errors.map((error, i) => (
            <div key={i} className="text-xs text-red-600">
              ✗ {error.message}
            </div>
          ))}
        </div>
      )}

      {hasWarnings && (
        <div className="space-y-1 mt-2">
          {validationResult.warnings.map((warning, i) => (
            <div key={i} className="text-xs text-yellow-600">
              ⚠ {warning.message}
            </div>
          ))}
        </div>
      )}
    </div>
  )
}

