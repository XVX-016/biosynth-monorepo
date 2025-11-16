import React, { useState, useMemo } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { getTemplates, loadTemplate, placeTemplate } from '../../kernel/templateLoader';
import { useTemplateToolStore } from '../../store/templateTool.store';
import { useUIStore } from '../../store/uiStore';
import { TemplateCategory } from './TemplateCategory';
import './TemplatePanel.css';

type TemplateCategoryType = 'rings' | 'groups' | 'atoms' | 'common';

interface TemplateWithCategory {
  id: string;
  name: string;
  category: TemplateCategoryType;
}

// Categorize templates
const categorizeTemplate = (templateId: string): TemplateCategoryType => {
  if (templateId === 'benzene') return 'rings';
  if (['water', 'ethanol'].includes(templateId)) return 'groups';
  if (templateId === 'methane') return 'atoms';
  return 'common';
};

// Simple fuzzy search
const fuzzyMatch = (text: string, query: string): boolean => {
  const textLower = text.toLowerCase();
  const queryLower = query.toLowerCase();
  
  // Exact match
  if (textLower.includes(queryLower)) return true;
  
  // Character sequence match (fuzzy)
  let textIndex = 0;
  for (let i = 0; i < queryLower.length; i++) {
    const char = queryLower[i];
    const foundIndex = textLower.indexOf(char, textIndex);
    if (foundIndex === -1) return false;
    textIndex = foundIndex + 1;
  }
  return true;
};

export const TemplatePanel: React.FC = () => {
  const [searchQuery, setSearchQuery] = useState('');
  const templates = getTemplates();
  const { startTemplateDrag, stopTemplateDrag } = useTemplateToolStore();
  const { templatePanelExpanded, toggleTemplatePanel } = useUIStore();

  // Add categories to templates
  const templatesWithCategories = useMemo<TemplateWithCategory[]>(() => {
    return templates.map((t) => ({
      id: t.id,
      name: t.name,
      category: categorizeTemplate(t.id),
    }));
  }, [templates]);

  // Filter templates based on search
  const filteredTemplates = useMemo(() => {
    if (!searchQuery.trim()) return templatesWithCategories;

    const query = searchQuery.trim();
    return templatesWithCategories.filter((t) =>
      fuzzyMatch(t.name, query) || fuzzyMatch(t.id, query)
    );
  }, [templatesWithCategories, searchQuery]);

  // Group by category
  const templatesByCategory = useMemo(() => {
    const grouped: Record<TemplateCategoryType, TemplateWithCategory[]> = {
      rings: [],
      groups: [],
      atoms: [],
      common: [],
    };

    filteredTemplates.forEach((template) => {
      grouped[template.category].push(template);
    });

    return grouped;
  }, [filteredTemplates]);

  const handleDragStart = (e: React.DragEvent, templateId: string) => {
    e.dataTransfer.effectAllowed = 'copy';
    e.dataTransfer.setData('template-id', templateId);
    startTemplateDrag(templateId);
  };

  const handleDragEnd = () => {
    stopTemplateDrag();
  };

  const handleTemplateClick = (templateId: string) => {
    const template = loadTemplate(templateId);
    if (template) {
      placeTemplate(template, { x: 0, y: 0, z: 0 });
    }
  };

  return (
    <div className="template-panel">
      {/* Header with collapse button */}
      <div className="template-panel-header">
        <h3 className="template-panel-title">Templates</h3>
        <button
          onClick={toggleTemplatePanel}
          className="template-panel-toggle"
          aria-label={templatePanelExpanded ? 'Collapse' : 'Expand'}
        >
          {templatePanelExpanded ? '▼' : '▶'}
        </button>
      </div>

      <AnimatePresence>
        {templatePanelExpanded && (
          <motion.div
            initial={{ opacity: 0, height: 0 }}
            animate={{ opacity: 1, height: 'auto' }}
            exit={{ opacity: 0, height: 0 }}
            className="template-panel-content"
          >
            {/* Search Bar */}
            <div className="template-panel-search">
              <input
                type="text"
                placeholder="Search templates..."
                value={searchQuery}
                onChange={(e) => setSearchQuery(e.target.value)}
                className="template-panel-search-input"
              />
            </div>

            {/* Template Categories */}
            <div className="template-panel-categories">
              <TemplateCategory
                title="Rings"
                templates={templatesByCategory.rings}
                onTemplateClick={handleTemplateClick}
                onDragStart={handleDragStart}
                onDragEnd={handleDragEnd}
              />
              <TemplateCategory
                title="Functional Groups"
                templates={templatesByCategory.groups}
                onTemplateClick={handleTemplateClick}
                onDragStart={handleDragStart}
                onDragEnd={handleDragEnd}
              />
              <TemplateCategory
                title="Core Atoms"
                templates={templatesByCategory.atoms}
                onTemplateClick={handleTemplateClick}
                onDragStart={handleDragStart}
                onDragEnd={handleDragEnd}
              />
              <TemplateCategory
                title="Common"
                templates={templatesByCategory.common}
                onTemplateClick={handleTemplateClick}
                onDragStart={handleDragStart}
                onDragEnd={handleDragEnd}
              />
            </div>

            {filteredTemplates.length === 0 && (
              <div className="template-panel-empty">
                No templates found
              </div>
            )}
          </motion.div>
        )}
      </AnimatePresence>
    </div>
  );
};
