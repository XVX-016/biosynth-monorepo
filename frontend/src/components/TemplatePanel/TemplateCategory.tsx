import React from 'react';
import { TemplateItem, TemplateItemProps } from './TemplateItem';

export interface TemplateCategoryProps {
  title: string;
  templates: Array<{
    id: string;
    name: string;
    category: string;
  }>;
  onTemplateClick: (templateId: string) => void;
  onDragStart: (e: React.DragEvent, templateId: string) => void;
  onDragEnd: () => void;
}

export const TemplateCategory: React.FC<TemplateCategoryProps> = ({
  title,
  templates,
  onTemplateClick,
  onDragStart,
  onDragEnd,
}) => {
  if (templates.length === 0) return null;

  return (
    <div className="mb-6">
      <h4 className="text-sm font-bold text-chrome mb-3 uppercase tracking-wider">
        {title}
      </h4>
      <div className="grid grid-cols-2 gap-2">
        {templates.map((template) => (
          <TemplateItem
            key={template.id}
            id={template.id}
            name={template.name}
            category={template.category}
            onDragStart={onDragStart}
            onDragEnd={onDragEnd}
            onClick={() => onTemplateClick(template.id)}
          />
        ))}
      </div>
    </div>
  );
};

