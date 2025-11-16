import { create } from 'zustand';

interface UIState {
  templatePanelExpanded: boolean;
  toggleTemplatePanel: () => void;
  setTemplatePanelExpanded: (expanded: boolean) => void;
}

export const useUIStore = create<UIState>((set) => ({
  templatePanelExpanded: true,
  toggleTemplatePanel: () => {
    set((state) => ({ templatePanelExpanded: !state.templatePanelExpanded }));
  },
  setTemplatePanelExpanded: (expanded: boolean) => {
    set({ templatePanelExpanded: expanded });
  },
}));

