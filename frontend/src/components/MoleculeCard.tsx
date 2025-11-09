import React from 'react'
import { motion } from 'framer-motion'
import { MoleculeItem } from '../lib/api'

interface MoleculeCardProps {
  item: MoleculeItem
  onOpen: () => void
  onDelete: () => void
}

export default function MoleculeCard({ item, onOpen, onDelete }: MoleculeCardProps) {
  const formatDate = (dateStr: string) => {
    const date = new Date(dateStr)
    return date.toLocaleDateString('en-US', { month: 'short', day: 'numeric', year: 'numeric' })
  }

  return (
    <motion.div
      initial={{ opacity: 0, y: 20 }}
      animate={{ opacity: 1, y: 0 }}
      className="bg-panel p-4 rounded-xl shadow-soft hover:shadow-lg transition-shadow"
    >
      {/* Thumbnail */}
      <div className="h-40 bg-aluminum-light rounded-lg flex items-center justify-center overflow-hidden mb-3">
        {item.thumbnail_b64 ? (
          <img
            src={item.thumbnail_b64}
            alt={item.name}
            className="max-h-full max-w-full object-contain"
          />
        ) : (
          <div className="text-text-tertiary text-sm">No preview</div>
        )}
      </div>

      {/* Info */}
      <h3 className="font-semibold text-text-primary mb-1 truncate">{item.name}</h3>
      {item.smiles && (
        <div className="text-xs text-text-secondary mb-2 font-mono truncate" title={item.smiles}>
          {item.smiles}
        </div>
      )}
      <div className="text-xs text-text-tertiary mb-3">{formatDate(item.created_at)}</div>

      {/* Actions */}
      <div className="flex gap-2">
        <motion.button
          whileHover={{ scale: 1.05 }}
          whileTap={{ scale: 0.95 }}
          onClick={onOpen}
          className="flex-1 px-3 py-2 rounded-lg bg-accent-blue text-white text-sm font-medium hover:opacity-90"
        >
          Open in Lab
        </motion.button>
        <motion.button
          whileHover={{ scale: 1.05 }}
          whileTap={{ scale: 0.95 }}
          onClick={onDelete}
          className="px-3 py-2 rounded-lg bg-aluminum-DEFAULT text-text-primary text-sm font-medium hover:bg-aluminum-dark"
        >
          Delete
        </motion.button>
      </div>
    </motion.div>
  )
}

