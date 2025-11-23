import React, { useState, Suspense, useEffect } from 'react'
import { motion } from 'framer-motion'
import { useInView } from 'react-intersection-observer'
import type { MoleculeItem } from '../lib/api'
import BarbellViewer from './BarbellViewer'
import { useGPUSafe } from '../hooks/useGPUSafe'
import { convertSMILESToMolfile, saveMolfile } from '../lib/api'
import { updateUserMolecule } from '../lib/userMoleculeStore'
import { supabase } from '../supabase'

interface MoleculeCardProps {
  item: MoleculeItem & { molfile?: string | null; formula?: string | null; user_id?: string }
  onOpen: () => void
  onDelete?: () => void
  onFork?: () => void // For public library molecules
  showFork?: boolean // Whether to show fork button
  onMolfileUpdated?: () => void // Callback when molfile is persisted
}

export default function MoleculeCard({ 
  item, 
  onOpen, 
  onDelete,
  onFork,
  showFork = false,
  onMolfileUpdated
}: MoleculeCardProps) {
  const [hovered, setHovered] = useState(false)
  const [convertedMolfile, setConvertedMolfile] = useState<string | null>(null)
  const [converting, setConverting] = useState(false)
  const [saving, setSaving] = useState(false)
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

  // Get molfile from item or converted molfile
  const molfile = item.molfile && item.molfile.trim().length > 0 
    ? item.molfile 
    : convertedMolfile

  const has3DData = !!(molfile && molfile.trim().length > 0)
  
  // Always show 3D when available and in view (not just on hover)
  const show3D = inView && has3DData && isGPUSafe && !converting
  
  // Show thumbnail as background only if no 3D data available
  const showThumbnail = !!(item.thumbnail_b64 && !has3DData && !converting)

  // Auto-convert SMILES to molfile when in view and molfile is missing
  useEffect(() => {
    if (inView && !molfile && item.smiles && item.smiles.trim() && !converting && !convertedMolfile) {
      setConverting(true)
      convertSMILESToMolfile(item.smiles)
        .then(async (result) => {
          if (result.molfile && result.molfile.trim().length > 0) {
            setConvertedMolfile(result.molfile)
            
            // Persist molfile to database
            try {
              setSaving(true)
              
              // Try to save via Supabase if user_id is available (user_molecules)
              if (item.user_id && typeof item.id === 'string') {
                await updateUserMolecule(item.user_id, item.id, { molfile: result.molfile })
              } 
              // Try backend API if numeric ID (legacy molecules table)
              else if (typeof item.id === 'number') {
                await saveMolfile(item.id, result.molfile)
              }
              // Try Supabase public_molecules if no user_id
              else if (supabase && typeof item.id === 'string') {
                const { error } = await supabase
                  .from('public_molecules')
                  .update({ molfile: result.molfile })
                  .eq('id', item.id)
                
                if (error) throw error
              }
              
              // Notify parent component if callback provided
              if (onMolfileUpdated) {
                onMolfileUpdated()
              }
            } catch (error) {
              console.warn('Failed to persist molfile:', error)
              // Continue anyway - molfile is cached in state
            } finally {
              setSaving(false)
            }
          }
        })
        .catch((error) => {
          console.warn('Failed to convert SMILES to molfile:', error)
          // Don't show error to user, just fall back to thumbnail
        })
        .finally(() => {
          setConverting(false)
        })
    }
  }, [inView, molfile, item.smiles, item.id, item.user_id, converting, convertedMolfile, onMolfileUpdated])
  
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
        
        {/* Loading skeleton while converting or saving */}
        {(converting || saving) && (
          <div className="absolute inset-0 z-10 flex items-center justify-center bg-gray-50/90 backdrop-blur-sm">
            <div className="flex flex-col items-center gap-2">
              <div className="w-8 h-8 border-4 border-blue-500 border-t-transparent rounded-full animate-spin"></div>
              <div className="text-xs text-gray-600">
                {converting ? 'Generating 3D preview...' : 'Saving...'}
              </div>
            </div>
          </div>
        )}

        {/* 3D Viewer - always visible when molfile is available, interactive on hover */}
        {show3D && molfile && !converting && (
          <div className="absolute inset-0 z-10">
            <Suspense fallback={
              <div className="w-full h-full flex items-center justify-center bg-gray-50">
                <div className="w-6 h-6 border-2 border-gray-400 border-t-transparent rounded-full animate-spin"></div>
              </div>
            }>
              <BarbellViewer
                molfile={molfile}
                mode="card"
                height={160}
                interactive={hovered}
                autorotate={hovered}
              />
            </Suspense>
          </div>
        )}
        
        {/* Fallback when no data available */}
        {!has3DData && !showThumbnail && !converting && (
          <div className="w-full h-full flex flex-col items-center justify-center text-midGrey text-sm bg-zinc-100 relative group">
            <div className="text-xs mb-1">No preview</div>
            <div className="text-xs opacity-60">{item.name}</div>
            <div className="absolute bottom-2 left-2 right-2 opacity-0 group-hover:opacity-100 transition-opacity">
              <div className="bg-black/80 text-white text-xs px-2 py-1 rounded text-center">
                {item.smiles 
                  ? '3D preview will be generated on hover if SMILES is available.'
                  : '3D preview requires molfile or SMILES data. Thumbnails are generated when molecules are saved.'}
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

