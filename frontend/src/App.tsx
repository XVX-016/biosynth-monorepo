import React, { useState } from 'react'
import Dashboard from './pages/Dashboard'
import Library from './pages/Library'

type Page = 'dashboard' | 'library'

export default function App(){
  const [currentPage, setCurrentPage] = useState<Page>('dashboard')

  return (
    <div className="min-h-screen bg-aluminum-light">
      {/* Navigation */}
      <nav className="bg-panel border-b border-aluminum-DEFAULT px-8 py-4">
        <div className="flex items-center justify-between">
          <div className="flex items-center gap-8">
            <h1 className="text-xl font-bold text-text-primary">BioSynth AI</h1>
            <div className="flex gap-4">
              <button
                onClick={() => setCurrentPage('dashboard')}
                className={`px-4 py-2 rounded-lg font-medium transition-colors ${
                  currentPage === 'dashboard'
                    ? 'bg-accent-blue text-white'
                    : 'text-text-secondary hover:text-text-primary'
                }`}
              >
                Lab
              </button>
              <button
                onClick={() => setCurrentPage('library')}
                className={`px-4 py-2 rounded-lg font-medium transition-colors ${
                  currentPage === 'library'
                    ? 'bg-accent-blue text-white'
                    : 'text-text-secondary hover:text-text-primary'
                }`}
              >
                Library
              </button>
            </div>
          </div>
        </div>
      </nav>

      {/* Page Content */}
      <div className="p-8">
        {currentPage === 'dashboard' && <Dashboard />}
        {currentPage === 'library' && <Library />}
      </div>
    </div>
  )
}
