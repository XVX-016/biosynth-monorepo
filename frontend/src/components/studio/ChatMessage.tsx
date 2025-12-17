import { motion } from 'framer-motion';
import type { StudioMessage } from '../../types/studio';

interface ChatMessageProps {
    message: StudioMessage;
}

export default function ChatMessage({ message }: ChatMessageProps) {
    const isMentor = message.role === 'mentor';

    return (
        <motion.div
            initial={{ opacity: 0, y: 10, scale: 0.95 }}
            animate={{ opacity: 1, y: 0, scale: 1 }}
            transition={{ duration: 0.4, ease: "easeOut" }}
            className={`flex ${isMentor ? 'justify-start' : 'justify-end'} mb-4`}
        >
            <div
                className={`max-w-[85%] px-4 py-3 rounded-2xl shadow-sm ${isMentor
                    ? 'bg-transparent text-black font-medium leading-relaxed' // Floating text for mentor
                    : 'bg-white/40 backdrop-blur-xl border border-white/20 text-black shadow-lg rounded-tr-none' // Glassmorphic for user
                    }`}
            >
                {isMentor && (
                    <div className="text-[10px] font-bold uppercase tracking-widest text-darkGrey/40 mb-1">
                        Mentor Response
                    </div>
                )}
                <div className="text-sm whitespace-pre-wrap">
                    {message.content}
                </div>
            </div>
        </motion.div>
    );
}
