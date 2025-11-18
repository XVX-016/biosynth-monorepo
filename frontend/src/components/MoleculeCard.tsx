import React from 'react'
import { motion } from 'framer-motion'
import type { MoleculeItem } from '../lib/api'

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
      className="bg-white p-4 rounded-xl shadow-neon border border-lightGrey hover:shadow-neon-hover hover:border-midGrey transition-all"
    >
      {/* Thumbnail */}
      <div className="h-40 bg-offwhite rounded-lg flex items-center justify-center overflow-hidden mb-3">
        {item.thumbnail_b64 ? (
          <img
            src={item.thumbnail_b64}
            alt={item.name}
            className="max-h-full max-w-full object-contain"
          />
        ) : (
          <div className="text-midGrey text-sm">No preview</div>
        )}
      </div>

      {/* Info */}
      <h3 className="font-semibold text-black mb-1 truncate">{item.name}</h3>
      {item.smiles && (
        <div className="text-xs text-darkGrey mb-2 font-mono truncate" title={item.smiles}>
          {item.smiles}
        </div>
      )}
      <div className="text-xs text-midGrey mb-3">{formatDate(item.created_at)}</div>

      {/* Actions */}
      <div className="flex gap-2">
        <motion.button
          whileHover={{ scale: 1.05 }}
          whileTap={{ scale: 0.95 }}
          onClick={onOpen}
          className="flex-1 px-3 py-2 rounded-lg bg-black text-white text-sm font-medium hover:bg-darkGrey hover:shadow-neon transition-all"
        >
          Open in Lab
        </motion.button>
        <motion.button
          whileHover={{ scale: 1.05 }}
          whileTap={{ scale: 0.95 }}
          onClick={onDelete}
          className="px-3 py-2 rounded-lg bg-white text-darkGrey hover:text-black hover:bg-offwhite border border-lightGrey text-sm font-medium transition-all"
        >
          Delete
        </motion.button>
      </div>
    </motion.div>
  )
}

