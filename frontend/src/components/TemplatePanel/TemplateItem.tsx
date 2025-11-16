import React from 'react';
import { motion } from 'framer-motion';

export interface TemplateItemProps {
  id: string;
  name: string;
  category: string;
  onDragStart: (e: React.DragEvent, templateId: string) => void;
  onDragEnd: () => void;
  onClick: () => void;
}

export const TemplateItem: React.FC<TemplateItemProps> = ({
  id,
  name,
  category,
  onDragStart,
  onDragEnd,
  onClick,
}) => {
  return (
    <motion.div
      draggable
      data-template-id={id}
      onDragStart={(e) => onDragStart(e, id)}
      onDragEnd={onDragEnd}
      whileHover={{ scale: 1.02 }}
      whileTap={{ scale: 0.98 }}
      className="cursor-grab active:cursor-grabbing"
    >
      <div
        onClick={onClick}
        className="p-3 rounded-lg bg-ionBlack/30 border border-chrome/30 hover:border-neonCyan/50 hover:shadow-neon-sm transition-all duration-200"
      >
        {/* Preview Icon - Circle Atom Graph */}
        <div className="w-full h-16 mb-2 flex items-center justify-center bg-frostedGlass rounded border border-chrome/20">
          <div className="relative w-12 h-12">
            {/* Simple atom representation */}
            <div className="absolute inset-0 flex items-center justify-center">
              <div className="w-2 h-2 rounded-full bg-neonCyan/60"></div>
            </div>
            <div className="absolute top-0 left-1/2 -translate-x-1/2 w-1 h-1 rounded-full bg-chrome/40"></div>
            <div className="absolute bottom-0 left-1/2 -translate-x-1/2 w-1 h-1 rounded-full bg-chrome/40"></div>
            <div className="absolute left-0 top-1/2 -translate-y-1/2 w-1 h-1 rounded-full bg-chrome/40"></div>
            <div className="absolute right-0 top-1/2 -translate-y-1/2 w-1 h-1 rounded-full bg-chrome/40"></div>
          </div>
        </div>

        {/* Template Name */}
        <div className="text-sm font-semibold text-ivory mb-1 truncate">{name}</div>

        {/* Category Tag */}
        <div className="text-xs text-spaceGrey bg-chrome/20 px-2 py-0.5 rounded inline-block">
          {category}
        </div>
      </div>
    </motion.div>
  );
};

