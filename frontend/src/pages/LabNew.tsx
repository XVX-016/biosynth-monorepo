import React from 'react'
import { motion } from 'framer-motion'
import LabViewer from '../components/lab/LabViewer'
import LabToolPanel from '../components/lab/LabToolPanel'

/**
 * New Lab page using the foundation architecture
 * This is a clean implementation with the new labStore
 */
export default function LabNew() {
  return (
    <motion.div 
      className="lab-layout p-4"
      initial={{ opacity: 0 }}
      animate={{ opacity: 1 }}
      transition={{ duration: 0.3 }}
    >
      <div className="grid grid-cols-12 gap-4 h-[calc(100vh-2rem)]">
        <motion.div 
          className="col-span-12 lg:col-span-3"
          initial={{ opacity: 0, x: -20 }}
          animate={{ opacity: 1, x: 0 }}
          transition={{ duration: 0.4, delay: 0.1 }}
        >
          <div className="bg-white rounded-xl shadow-lg border border-gray-200 p-4 h-full overflow-y-auto">
            <LabToolPanel />
          </div>
        </motion.div>
        <motion.div 
          className="col-span-12 lg:col-span-9"
          initial={{ opacity: 0, scale: 0.98 }}
          animate={{ opacity: 1, scale: 1 }}
          transition={{ duration: 0.4, delay: 0.2 }}
        >
          <div className="relative rounded-xl shadow-lg border border-gray-200 p-2 h-full bg-gray-50">
            <LabViewer />
          </div>
        </motion.div>
      </div>
    </motion.div>
  )
}

