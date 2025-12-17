import { useState } from 'react';
import { motion } from 'framer-motion';
import MentorStrip from '../components/studio/MentorStrip';
import ChatInterface from '../components/studio/ChatInterface';
import Studio3DScene from '../components/studio/Studio3DScene';
import Card from '../components/ui/Card';
import type { StudioMode } from '../types/studio';
import { MENTORS } from '../config/mentors';
import type { MoleculeGraph } from '@biosynth/engine';

export default function StudioPage() {
  const [activeMentorId, setActiveMentorId] = useState<string | undefined>(MENTORS[0].id);
  const [currentMode, setCurrentMode] = useState<StudioMode>('design');
  const [currentMolecule] = useState<MoleculeGraph | null>(null);

  const activeMentor = MENTORS.find(m => m.id === activeMentorId);

  return (
    <motion.div
      initial={{ opacity: 0 }}
      animate={{ opacity: 1 }}
      transition={{ duration: 0.3 }}
      className="space-y-8 max-w-7xl mx-auto"
    >
      {/* Studio Header */}
      <div className="space-y-2">
        <h1 className="text-4xl font-extrabold text-black tracking-tight">MolForge Studio</h1>
        <p className="text-darkGrey text-lg">
          Design, optimize, and simulate molecules with AI chemistry mentors
        </p>
      </div>

      {/* Mentor Selection Strip */}
      <MentorStrip
        activeMentorId={activeMentorId}
        onMentorSelect={setActiveMentorId}
      />

      {/* Main Workspace */}
      <div className="grid grid-cols-1 lg:grid-cols-12 gap-6 min-h-[600px]">
        {/* Intent Panel (Left) */}
        <div className="lg:col-span-4 space-y-4 flex flex-col">
          <Card className="flex-1 p-6 bg-white flex flex-col">
            <div className="flex items-center justify-between mb-6">
              <h2 className="text-xl font-bold text-black">Intent Panel</h2>
              {/* Mode Selector */}
              <div className="flex bg-offwhite p-1 rounded-lg border border-lightGrey">
                {(['design', 'optimize', 'simulate'] as StudioMode[]).map((mode) => (
                  <button
                    key={mode}
                    onClick={() => setCurrentMode(mode)}
                    className={`px-3 py-1.5 text-xs font-semibold rounded-md transition-all ${currentMode === mode
                      ? 'bg-black text-white shadow-sm'
                      : 'text-darkGrey hover:text-black'
                      }`}
                  >
                    {mode.charAt(0).toUpperCase() + mode.slice(1)}
                  </button>
                ))}
              </div>
            </div>

            {/* Chat/Input Area */}
            <div className="flex-1 flex flex-col overflow-hidden">
              <ChatInterface
                mentorId={activeMentorId!}
                mode={currentMode}
              />
            </div>
          </Card>
        </div>

        {/* Molecule Canvas (Right) */}
        <div className="lg:col-span-8 relative">
          <Card className="h-full overflow-hidden bg-offwhite border-none relative flex items-center justify-center">
            {/* Background Mentor Visual */}
            <div className="absolute inset-0">
              <Studio3DScene
                mentorId={activeMentorId!}
                accentColor={activeMentor?.accentColor || '#000000'}
                modelPath={activeMentor?.avatarModelPath || ''}
                molecule={currentMolecule}
              />
            </div>

            {/* UI Overlays */}
            <div className="absolute top-6 left-6 z-10">
              <div className="bg-white/80 backdrop-blur-md border border-lightGrey px-4 py-2 rounded-full flex items-center gap-2 shadow-sm">
                <div className="w-2 h-2 rounded-full animate-pulse" style={{ background: activeMentor?.accentColor }} />
                <span className="text-xs font-bold text-black uppercase tracking-wider">{activeMentor?.name}</span>
                <span className="text-xs text-darkGrey">• Active</span>
              </div>
            </div>

            <div className="z-10 text-center">
              <div className="text-darkGrey/40 font-mono text-sm mb-4">MOLECULE_CANVAS_V1</div>
              <div className="w-64 h-64 border-2 border-dashed border-lightGrey rounded-full flex items-center justify-center">
                <div className="text-midGrey text-xs text-center p-8">
                  3D Visualization and Liquid Chat overlays will be integrated here.
                </div>
              </div>
            </div>
          </Card>
        </div>
      </div>
    </motion.div>
  );
}
