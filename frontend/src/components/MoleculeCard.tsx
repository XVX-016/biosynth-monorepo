import { useState, Suspense } from 'react'
import { motion } from 'framer-motion'
import { useInView } from 'react-intersection-observer'
import type { MoleculeItem } from '../lib/api'
import BarbellViewer from './BarbellViewer'
import { useGPUSafe } from '../hooks/useGPUSafe'

interface MoleculeCardProps {
  item: MoleculeItem & { molfile?: string | null; formula?: string | null; user_id?: string }
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
  
  // Lazy load 3D viewer when card enters viewport
  const { ref, inView } = useInView({
    triggerOnce: true,
    threshold: 0.1,
    rootMargin: '50px',
  })

  const formatDate = (dateStr: string) => {
    const date = new Date(dateStr)
    return date.toLocaleDateString('en-US', { month: 'short', day: 'numeric', year: 'numeric' })
  }

  // Only use molfile from item - NO conversion in card
  const molfile = item.molfile && item.molfile.trim().length > 0 ? item.molfile : null
  const has3DData = !!molfile
  
  // Always show 3D when available and in view
  const show3D = inView && has3DData && isGPUSafe
  
  // Show thumbnail as background only if no 3D data available
  const showThumbnail = !!(item.thumbnail_b64 && !has3DData)

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
      {/* 3D Viewer or Thumbnail with lazy loading - FIXED HEIGHT to prevent re-renders */}
      <div className="h-[180px] bg-offwhite rounded-lg overflow-hidden mb-3 relative">
        {/* Background thumbnail (shown only as fallback when no 3D data) */}
        {showThumbnail && item.thumbnail_b64 && (
          <img
            src={item.thumbnail_b64.startsWith('data:') 
              ? item.thumbnail_b64 
              : `data:image/png;base64,${item.thumbnail_b64}`}
            alt={item.name}
            className="w-full h-full object-contain"
            onError={(e) => {
              // Hide broken image
              (e.target as HTMLImageElement).style.display = 'none';
            }}
          />
        )}
        
        {/* 3D Viewer - always visible when molfile is available, interactive on hover */}
        {show3D && molfile && (
          <div className="absolute inset-0 z-10" key={`mol-${item.id}`}>
            <Suspense fallback={
              <div className="w-full h-full flex items-center justify-center bg-gray-50">
                <div className="w-6 h-6 border-2 border-gray-400 border-t-transparent rounded-full animate-spin"></div>
              </div>
            }>
              <BarbellViewer
                molfile={molfile}
                mode="card"
                height={180}
                hovered={hovered}
                interactive={false}
                autorotate={hovered}
                atomScale={0.28}
                bondRadius={0.06}
              />
            </Suspense>
          </div>
        )}
        
        {/* Fallback when no data available */}
        {!has3DData && !showThumbnail && (
          <div className="w-full h-full flex flex-col items-center justify-center text-midGrey text-sm bg-zinc-100 relative group">
            <div className="text-xs mb-1">No preview</div>
            <div className="text-xs opacity-60">{item.name}</div>
            <div className="absolute bottom-2 left-2 right-2 opacity-0 group-hover:opacity-100 transition-opacity">
              <div className="bg-black/80 text-white text-xs px-2 py-1 rounded text-center">
                3D preview requires molfile data. Conversion happens automatically when molecules are loaded.
              </div>
            </div>
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

