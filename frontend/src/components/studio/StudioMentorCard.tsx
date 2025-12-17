import { motion } from 'framer-motion';
import type { StudioMentor } from '../../types/studio';
import Card from '../ui/Card';

interface StudioMentorCardProps {
    mentor: StudioMentor;
    isActive?: boolean;
    onClick: (mentorId: string) => void;
}

export default function StudioMentorCard({ mentor, isActive, onClick }: StudioMentorCardProps) {
    return (
        <motion.div
            whileHover={{ scale: 1.02, y: -5 }}
            whileTap={{ scale: 0.98 }}
            className="relative cursor-pointer"
            onClick={() => onClick(mentor.id)}
        >
            <Card
                className={`w-full overflow-hidden transition-all duration-300 border-2 ${isActive
                    ? 'border-black shadow-neon ring-1 ring-black/10'
                    : 'border-transparent hover:border-black/20'
                    }`}
                style={{ background: mentor.backgroundGradient }}
            >
                <div className="flex flex-col h-48 relative p-6">
                    {/* Header Info */}
                    <div className="z-10">
                        <h3 className="text-xl font-bold text-white mb-0.5">{mentor.name}</h3>
                        <p className="text-white/80 text-sm font-medium">{mentor.role}</p>
                    </div>

                    {/* Upper Body Avatar Placeholder (3D Canvas will go here later) */}
                    <div className="absolute inset-0 flex items-end justify-center pointer-events-none">
                        <div className="w-full h-full bg-gradient-to-t from-black/20 to-transparent absolute bottom-0" />
                        {/* 3D Model will be rendered by a dedicated component hooked into this layout */}
                        <div className="w-40 h-40 bg-white/10 rounded-full blur-xl mb-[-20%] opacity-50" />
                    </div>

                    {/* Try Mentor Hook */}
                    <div className="absolute bottom-6 right-6 z-10">
                        <span className="text-white/90 text-sm font-semibold flex items-center gap-1">
                            Try Mentor <span className="text-lg">→</span>
                        </span>
                    </div>
                </div>
            </Card>

            {/* Active Glow Ring */}
            {isActive && (
                <motion.div
                    layoutId="activeGlow"
                    className="absolute inset-0 border-2 border-black rounded-xl -m-1 pointer-events-none"
                    initial={false}
                    transition={{ type: "spring", stiffness: 300, damping: 30 }}
                />
            )}
        </motion.div>
    );
}
