import React, { useState, useMemo } from 'react';
import { motion } from 'framer-motion';
import { getTemplates } from '../../kernel/templateLoader';
import { useTemplateToolStore } from '../../store/templateTool.store';
import Button from '../ui/Button';

type TemplateCategory = 'all' | 'rings' | 'groups' | 'atoms';

export const TemplatePanel: React.FC = () => {
  const [searchQuery, setSearchQuery] = useState('');
  const [activeCategory, setActiveCategory] = useState<TemplateCategory>('all');
  const templates = getTemplates();
  const { startTemplateDrag, stopTemplateDrag } = useTemplateToolStore();

  // Filter templates based on search and category
  const filteredTemplates = useMemo(() => {
    let filtered = templates;

    // Category filter
    if (activeCategory === 'rings') {
      filtered = filtered.filter((t) => t.id === 'benzene');
    } else if (activeCategory === 'groups') {
      filtered = filtered.filter((t) => ['water', 'ethanol'].includes(t.id));
    } else if (activeCategory === 'atoms') {
      filtered = filtered.filter((t) => t.id === 'methane');
    }

    // Search filter
    if (searchQuery.trim()) {
      const query = searchQuery.toLowerCase();
      filtered = filtered.filter(
        (t) =>
          t.name.toLowerCase().includes(query) ||
          t.id.toLowerCase().includes(query)
      );
    }

    return filtered;
  }, [templates, searchQuery, activeCategory]);

  const handleDragStart = (e: React.DragEvent, templateId: string) => {
    e.dataTransfer.effectAllowed = 'copy';
    e.dataTransfer.setData('template-id', templateId);
    startTemplateDrag(templateId);
  };

  const handleDragEnd = () => {
    stopTemplateDrag();
  };

  const categories: Array<{ id: TemplateCategory; label: string }> = [
    { id: 'all', label: 'All' },
    { id: 'rings', label: 'Rings' },
    { id: 'groups', label: 'Functional Groups' },
    { id: 'atoms', label: 'Core Atoms' },
  ];

  return (
    <div className="p-4 bg-frostedGlass rounded-lg border border-chrome/30 w-full">
      <h3 className="text-lg font-bold mb-3 text-ivory">Templates</h3>

      {/* Search Bar */}
      <div className="mb-4">
        <input
          type="text"
          placeholder="Search templates..."
          value={searchQuery}
          onChange={(e) => setSearchQuery(e.target.value)}
          className="w-full px-3 py-2 rounded-lg bg-ionBlack/50 border border-chrome/30 text-ivory placeholder-spaceGrey focus:outline-none focus:border-neonCyan/50 focus:ring-2 focus:ring-neonCyan/20"
        />
      </div>

      {/* Category Tabs */}
      <div className="flex gap-2 mb-4 flex-wrap">
        {categories.map((cat) => (
          <button
            key={cat.id}
            onClick={() => setActiveCategory(cat.id)}
            className={`px-3 py-1 rounded text-sm font-medium transition-all ${
              activeCategory === cat.id
                ? 'bg-chrome text-spaceGrey border border-neonCyan/50'
                : 'bg-ionBlack/30 text-ivory border border-chrome/20 hover:border-chrome/40'
            }`}
          >
            {cat.label}
          </button>
        ))}
      </div>

      {/* Template Grid */}
      <div className="grid grid-cols-2 gap-2 max-h-96 overflow-y-auto">
        {filteredTemplates.length === 0 ? (
          <div className="col-span-2 text-center text-spaceGrey py-4">
            No templates found
          </div>
        ) : (
          filteredTemplates.map((template) => (
            <motion.div
              key={template.id}
              draggable
              data-template-id={template.id}
              onDragStart={(e) => handleDragStart(e as any, template.id)}
              onDragEnd={handleDragEnd}
              whileHover={{ scale: 1.02 }}
              whileTap={{ scale: 0.98 }}
              className="cursor-grab active:cursor-grabbing"
            >
              <Button
                variant="secondary"
                className="w-full text-left justify-start"
                onClick={() => {
                  // Click handler can also trigger template placement
                  // This will be handled by the drag-drop system
                }}
              >
                <div className="flex flex-col items-start">
                  <span className="font-semibold text-ivory">{template.name}</span>
                  <span className="text-xs text-spaceGrey mt-1">
                    {template.id}
                  </span>
                </div>
              </Button>
            </motion.div>
          ))
        )}
      </div>
    </div>
  );
};

