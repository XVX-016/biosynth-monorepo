import React, { useState } from 'react'
import { motion } from 'framer-motion'
import Navbar from '../../components/Navbar'
import LabToolPanel from '../../components/lab/LabToolPanel'
import LabViewer from '../../components/lab/LabViewer'

/**
 * Clean white minimal Lab layout with standard navigation
 */
export default function NewLabLayout() {
  const [mobileOpen, setMobileOpen] = useState(false)

  return (
    <div className="h-screen w-screen flex flex-col bg-white overflow-hidden">
      {/* Standard Navbar */}
      <Navbar onToggleMenu={() => setMobileOpen((v) => !v)} />

      {/* Main Content Area */}
      <div className="flex flex-1 overflow-hidden">
        {/* Left Sidebar - Tools */}
        <motion.aside
          initial={{ x: -20, opacity: 0 }}
          animate={{ x: 0, opacity: 1 }}
          transition={{ duration: 0.4, ease: 'easeOut' }}
          className="w-20 border-r border-[#dcdcdc] bg-[#f8f9fa] flex flex-col shrink-0 overflow-y-auto"
        >
          <LabToolPanel />
        </motion.aside>

        {/* Main 3D Workspace */}
        <motion.main
          initial={{ opacity: 0, scale: 0.98 }}
          animate={{ opacity: 1, scale: 1 }}
          transition={{ duration: 0.5 }}
          className="flex-1 bg-white relative overflow-hidden"
        >
          {/* Soft vignette effect */}
          <div className="pointer-events-none absolute inset-0 bg-gradient-radial from-transparent via-transparent to-black/5 z-10" />
          <LabViewer />
        </motion.main>
      </div>

      {/* Bottom Status Bar */}
      <motion.footer
        initial={{ y: 20, opacity: 0 }}
        animate={{ y: 0, opacity: 1 }}
        transition={{ duration: 0.4, delay: 0.2 }}
        className="h-8 border-t border-[#dcdcdc] text-xs text-gray-500 flex items-center justify-between px-4 bg-white shrink-0"
      >
        <div className="flex items-center gap-4">
          <span>Ready</span>
        </div>
        <div className="text-[10px] text-gray-400">
          Press spacebar for help
        </div>
      </motion.footer>
    </div>
  )
}
