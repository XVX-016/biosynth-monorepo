import React, { useState, Suspense } from 'react'
import { motion } from 'framer-motion'
import { useInView } from 'react-intersection-observer'
import type { MoleculeItem } from '../lib/api'
import BarbellViewer from './BarbellViewer'
import { useGPUSafe } from '../hooks/useGPUSafe'

interface MoleculeCardProps {
  item: MoleculeItem & { molfile?: string | null; formula?: string | null }
  onOpen: () => void
  onDelete?: () => void
  onFork?: () => void // For public library molecules
  showFork?: boolean // Whether to show fork button
}

export default function MoleculeCard({ 
  item, 
  onOpen, 
  onDelete,
  onFork,
  showFork = false
}: MoleculeCardProps) {
  const [hovered, setHovered] = useState(false)
  const isGPUSafe = useGPUSafe()
  
  // Lazy load 3D viewer only when card enters viewport and is hovered
  const { ref, inView } = useInView({
    triggerOnce: true,
    threshold: 0.1,
    rootMargin: '50px',
  })

  const formatDate = (dateStr: string) => {
    const date = new Date(dateStr)
    return date.toLocaleDateString('en-US', { month: 'short', day: 'numeric', year: 'numeric' })
  }

  const has3DData = !!(item.molfile && item.molfile.trim().length > 0)
  // Only show 3D on desktop devices (GPU safe) and when hovered/in view
  const show3D = hovered && inView && has3DData && isGPUSafe
  // Show thumbnail as background if it exists, or as fallback if no 3D data
  const showThumbnail = !!(item.thumbnail_b64 && (!has3DData || !show3D))
  
  // Debug logging (only log once per molecule to reduce noise)
  React.useEffect(() => {
    if (!item.molfile && !item.thumbnail_b64) {
      console.log('MoleculeCard - No preview data:', {
        name: item.name,
        molfile: item.molfile,
        thumbnail: !!item.thumbnail_b64,
        has3DData
      });
    }
  }, [item.name]); // Only log when name changes

  return (
    <motion.div
      ref={ref}
      layout
      initial={{ opacity: 0, y: 20 }}
      animate={{ opacity: 1, y: 0 }}
      whileHover={{ y: -6, boxShadow: '0px 10px 30px rgba(0,0,0,0.10)' }}
      onMouseEnter={() => setTimeout(() => setHovered(true), 150)}
      onMouseLeave={() => setHovered(false)}
      className="bg-white p-4 rounded-xl shadow-neon border border-lightGrey hover:shadow-neon-hover hover:border-midGrey transition-all cursor-pointer"
      onClick={onOpen}
    >
      {/* 3D Viewer or Thumbnail with lazy loading */}
      <div className="h-40 bg-offwhite rounded-lg overflow-hidden mb-3 relative">
        {/* Background thumbnail (shown when not hovering or as fallback) */}
        {showThumbnail && item.thumbnail_b64 && (
          <img
            src={item.thumbnail_b64.startsWith('data:') 
              ? item.thumbnail_b64 
              : `data:image/png;base64,${item.thumbnail_b64}`}
            alt={item.name}
            className={`w-full h-full object-contain transition-opacity duration-300 ${
              show3D ? 'opacity-20 blur-sm' : 'opacity-100'
            }`}
            onError={(e) => {
              // Hide broken image
              (e.target as HTMLImageElement).style.display = 'none';
            }}
          />
        )}
        
        {/* Lazy-mounted 3D Viewer - only renders when hovered and in viewport */}
        {show3D && item.molfile && (
          <motion.div
            className="absolute inset-0 z-10"
            initial={{ opacity: 0 }}
            animate={{ opacity: 1 }}
            transition={{ duration: 0.3 }}
            key={`3d-${item.id}-${hovered}`}
          >
            <Suspense fallback={
              <div className="w-full h-full flex items-center justify-center bg-gray-50">
                <div className="text-xs text-gray-400">Loading 3D...</div>
              </div>
            }>
              <BarbellViewer
                molfile={item.molfile}
                mode="card"
                height={160}
              />
            </Suspense>
          </motion.div>
        )}
        
        {/* Fallback when no data available */}
        {!has3DData && !showThumbnail && (
          <div className="w-full h-full flex flex-col items-center justify-center text-midGrey text-sm bg-zinc-100">
            <div className="text-xs mb-1">No preview</div>
            <div className="text-xs opacity-60">{item.name}</div>
          </div>
        )}
      </div>

      {/* Info */}
      <h3 className="font-semibold text-lg text-black mb-1 truncate">{item.name}</h3>
      {item.formula && (
        <div className="text-sm text-darkGrey mb-1">
          Formula: <span className="font-mono">{item.formula}</span>
        </div>
      )}
      {item.smiles && (
        <div className="text-xs text-midGrey mb-2 font-mono truncate" title={item.smiles}>
          SMILES: {item.smiles}
        </div>
      )}
      <div className="text-xs text-midGrey mb-3">{formatDate(item.created_at)}</div>

      {/* Actions */}
      <div className="flex gap-2" onClick={(e) => e.stopPropagation()}>
        {showFork && onFork && (
          <motion.button
            whileHover={{ scale: 1.05 }}
            whileTap={{ scale: 0.95 }}
            onClick={onFork}
            className="px-3 py-2 rounded-lg bg-blue-600 text-white text-sm font-medium hover:bg-blue-700 transition-all"
          >
            Fork
          </motion.button>
        )}
        <motion.button
          whileHover={{ scale: 1.05 }}
          whileTap={{ scale: 0.95 }}
          onClick={onOpen}
          className={`${showFork ? 'flex-1' : ''} px-3 py-2 rounded-lg bg-black text-white text-sm font-medium hover:bg-darkGrey hover:shadow-neon transition-all`}
        >
          Open in Lab
        </motion.button>
        {onDelete && (
          <motion.button
            whileHover={{ scale: 1.05 }}
            whileTap={{ scale: 0.95 }}
            onClick={onDelete}
            className="px-3 py-2 rounded-lg bg-white text-darkGrey hover:text-black hover:bg-offwhite border border-lightGrey text-sm font-medium transition-all"
          >
            Delete
          </motion.button>
        )}
      </div>
    </motion.div>
  )
}

