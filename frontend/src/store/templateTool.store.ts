import { create } from 'zustand';

interface TemplateToolState {
  isDraggingTemplate: boolean;
  draggedTemplateId: string | null;
  startTemplateDrag: (templateId: string) => void;
  stopTemplateDrag: () => void;
}

export const useTemplateToolStore = create<TemplateToolState>((set) => ({
  isDraggingTemplate: false,
  draggedTemplateId: null,
  startTemplateDrag: (templateId: string) => {
    set({ isDraggingTemplate: true, draggedTemplateId: templateId });
  },
  stopTemplateDrag: () => {
    set({ isDraggingTemplate: false, draggedTemplateId: null });
  },
}));

