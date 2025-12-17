import React, { useEffect, useRef } from 'react';
import { motion } from 'framer-motion';
import BarbellViewer from '../components/BarbellViewer';

// Caffeine MOLFILE
const CAFFEINE_MOLFILE = `
  MJ211200                      

 14 15  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  2  0  0  0  0
  3  4  1  0  0  0  0
  4  5  2  0  0  0  0
  5  6  1  0  0  0  0
  6  1  2  0  0  0  0
  1  7  1  0  0  0  0
  3  8  1  0  0  0  0
  5  9  1  0  0  0  0
  2 11  2  0  0  0  0
  4 12  2  0  0  0  0
  7 13  1  0  0  0  0
  8 14  1  0  0  0  0
  9 10  1  0  0  0  0
 10  6  1  0  0  0  0
M  END
`;

export default function Welcome() {
  const sceneRef = useRef<HTMLDivElement>(null);
  const [showScene, setShowScene] = React.useState(false);

  useEffect(() => {
    const handleScroll = () => {
      if (sceneRef.current) {
        const rect = sceneRef.current.getBoundingClientRect();
        const isVisible = rect.top < window.innerHeight && rect.bottom > 0;
        setShowScene(isVisible);
      }
    };

    // Check on mount for mobile
    if (window.innerWidth < 640) {
      handleScroll();
    } else {
      setShowScene(true);
    }

    window.addEventListener('scroll', handleScroll);
    return () => window.removeEventListener('scroll', handleScroll);
  }, []);


  return (
    <div className="min-h-screen bg-ionBlack text-ivory">
      <div className="grid grid-cols-1 lg:grid-cols-2 gap-8 p-6 lg:p-12">
        {/* Left Section */}
        <motion.div
          initial={{ opacity: 0, x: -50 }}
          animate={{ opacity: 1, x: 0 }}
          transition={{ duration: 0.6, ease: 'easeOut' }}
          className="flex flex-col justify-center space-y-8"
        >
          <motion.h1
            initial={{ opacity: 0, y: -20 }}
            animate={{ opacity: 1, y: 0 }}
            transition={{ duration: 0.6, delay: 0.2 }}
            className="text-5xl lg:text-6xl font-bold text-ivory"
          >
            Welcome to SynthCore
          </motion.h1>

          <motion.p
            initial={{ opacity: 0, y: -20 }}
            animate={{ opacity: 1, y: 0 }}
            transition={{ duration: 0.6, delay: 0.4 }}
            className="text-xl text-chrome leading-relaxed"
          >
            Design, simulate, and explore molecular structures with cutting-edge AI-powered tools.
            Build your molecular library and unlock the future of synthetic chemistry.
          </motion.p>

          <motion.div
            initial={{ opacity: 0, y: -20 }}
            animate={{ opacity: 1, y: 0 }}
            transition={{ duration: 0.6, delay: 0.6 }}
            className="flex flex-col sm:flex-row gap-4"
          >
            <motion.button
              whileHover={{ scale: 1.05 }}
              whileTap={{ scale: 0.95 }}
              className="px-8 py-4 rounded-lg font-semibold text-spaceGrey transition-all duration-300"
              style={{
                background: 'linear-gradient(135deg, #E3E6EB, #C0C5D2)',
                border: '2px solid',
                borderColor: 'rgba(139, 243, 255, 0.5)',
                boxShadow: '0 0 20px rgba(139, 243, 255, 0.3)',
              }}
              onHoverStart={(e) => {
                e.currentTarget.style.boxShadow = '0 0 30px rgba(139, 243, 255, 0.5)';
                e.currentTarget.style.borderColor = 'rgba(139, 243, 255, 0.8)';
              }}
              onHoverEnd={(e) => {
                e.currentTarget.style.boxShadow = '0 0 20px rgba(139, 243, 255, 0.3)';
                e.currentTarget.style.borderColor = 'rgba(139, 243, 255, 0.5)';
              }}
            >
              Get Started
            </motion.button>

            <motion.button
              whileHover={{ scale: 1.05 }}
              whileTap={{ scale: 0.95 }}
              className="px-8 py-4 rounded-lg font-semibold text-spaceGrey transition-all duration-300"
              style={{
                background: 'linear-gradient(135deg, #E3E6EB, #C0C5D2)',
                border: '2px solid',
                borderColor: 'rgba(139, 243, 255, 0.5)',
                boxShadow: '0 0 20px rgba(139, 243, 255, 0.3)',
              }}
              onHoverStart={(e) => {
                e.currentTarget.style.boxShadow = '0 0 30px rgba(139, 243, 255, 0.5)';
                e.currentTarget.style.borderColor = 'rgba(139, 243, 255, 0.8)';
              }}
              onHoverEnd={(e) => {
                e.currentTarget.style.boxShadow = '0 0 20px rgba(139, 243, 255, 0.3)';
                e.currentTarget.style.borderColor = 'rgba(139, 243, 255, 0.5)';
              }}
            >
              Learn More
            </motion.button>
          </motion.div>
        </motion.div>

        {/* Right Section */}
        <motion.div
          ref={sceneRef}
          initial={{ opacity: 0, x: 50 }}
          animate={{ opacity: showScene ? 1 : 0, x: showScene ? 0 : 50 }}
          transition={{ duration: 0.8, ease: 'easeOut' }}
          className="relative h-[600px] lg:h-full min-h-[600px]"
        >
          <div
            className="w-full h-full rounded-xl relative overflow-hidden"
            style={{
              background: 'rgba(255, 255, 255, 0.06)',
              backdropFilter: 'blur(10px)',
              border: '2px solid',
              borderImage: 'linear-gradient(135deg, #E3E6EB, #C0C5D2) 1',
              boxShadow: '0 8px 32px 0 rgba(15, 17, 21, 0.37)',
            }}
          >
            {/* FrostedGlass overlay with plasmaNeon edge */}
            <div
              className="absolute inset-0 pointer-events-none"
              style={{
                background: 'rgba(255, 255, 255, 0.06)',
                backdropFilter: 'blur(10px)',
                border: '1px solid rgba(192, 197, 210, 0.2)',
              }}
            />
            <div
              className="absolute inset-0 pointer-events-none"
              style={{
                background: 'linear-gradient(135deg, rgba(59, 199, 201, 0.1), rgba(139, 243, 255, 0.05))',
                boxShadow: 'inset 0 0 40px rgba(59, 199, 201, 0.2)',
              }}
            />

            {/* Hero Molecule 3D Viewer */}
            <div className="relative w-full h-full z-10 flex items-center justify-center">
              <BarbellViewer
                molfile={CAFFEINE_MOLFILE}
                mode="hero"
                height={600}
                className="w-full h-full bg-transparent"
              />
            </div>
          </div>
        </motion.div>
      </div>
    </div>
  );
}


