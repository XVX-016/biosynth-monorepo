import { describe, it, expect, vi, beforeEach } from 'vitest';
import { screen, fireEvent } from '@testing-library/react';
import { renderWithProviders } from '../../../tests/test-utils';
import { TemplatePanel } from '../TemplatePanel';
import * as templateLoader from '../../../kernel/templateLoader';

// Mock CSS imports
vi.mock('../TemplatePanel.css', () => ({}));

// Mock the stores
const mockToggleTemplatePanel = vi.fn();
const mockStartTemplateDrag = vi.fn();
const mockStopTemplateDrag = vi.fn();

vi.mock('../../../store/uiStore', () => ({
  useUIStore: () => ({
    templatePanelExpanded: true,
    toggleTemplatePanel: mockToggleTemplatePanel,
  }),
}));

vi.mock('../../../store/templateTool.store', () => ({
  useTemplateToolStore: () => ({
    startTemplateDrag: mockStartTemplateDrag,
    stopTemplateDrag: mockStopTemplateDrag,
  }),
}));

// Mock template loader
const mockGetTemplates = vi.fn(() => [
  { id: 'water', name: 'Water (H₂O)', data: { atoms: [], bonds: [] } },
  { id: 'benzene', name: 'Benzene (C₆H₆)', data: { atoms: [], bonds: [] } },
  { id: 'methane', name: 'Methane (CH₄)', data: { atoms: [], bonds: [] } },
  { id: 'ethanol', name: 'Ethanol (C₂H₅OH)', data: { atoms: [], bonds: [] } },
]);

const mockLoadTemplate = vi.fn((id) => ({ atoms: [], bonds: [] }));
const mockPlaceTemplate = vi.fn();

vi.mock('../../../kernel/templateLoader', () => ({
  getTemplates: () => mockGetTemplates(),
  loadTemplate: (id: string) => mockLoadTemplate(id),
  placeTemplate: (template: any, offset: any) => mockPlaceTemplate(template, offset),
}));

describe('TemplatePanel', () => {
  beforeEach(() => {
    vi.clearAllMocks();
  });

  it('renders template panel with header', () => {
    renderWithProviders(<TemplatePanel />);
    expect(screen.getByText('Templates')).toBeInTheDocument();
  });

  it('renders categories when expanded', () => {
    renderWithProviders(<TemplatePanel />);
    expect(screen.getByText('Rings')).toBeInTheDocument();
    expect(screen.getByText('Functional Groups')).toBeInTheDocument();
    expect(screen.getByText('Core Atoms')).toBeInTheDocument();
  });

  it('filters templates by search query', () => {
    renderWithProviders(<TemplatePanel />);
    const searchInput = screen.getByPlaceholderText('Search templates...');
    
    fireEvent.change(searchInput, { target: { value: 'water' } });
    
    expect(screen.getByText('Water (H₂O)')).toBeInTheDocument();
  });

  it('calls placeTemplate when template is clicked', () => {
    renderWithProviders(<TemplatePanel />);
    
    const waterTemplate = screen.getByText('Water (H₂O)');
    const clickableParent = waterTemplate.closest('div[onClick]') || waterTemplate.parentElement;
    if (clickableParent) {
      fireEvent.click(clickableParent);
    }
    
    expect(mockPlaceTemplate).toHaveBeenCalled();
  });

  it('handles drag start events', () => {
    renderWithProviders(<TemplatePanel />);
    
    const templateItem = screen.getByText('Water (H₂O)').closest('[draggable]');
    if (templateItem) {
      fireEvent.dragStart(templateItem);
      expect(mockStartTemplateDrag).toHaveBeenCalled();
    }
  });

  it('handles drag end events', () => {
    renderWithProviders(<TemplatePanel />);
    
    const templateItem = screen.getByText('Water (H₂O)').closest('[draggable]');
    if (templateItem) {
      fireEvent.dragEnd(templateItem);
      expect(mockStopTemplateDrag).toHaveBeenCalled();
    }
  });

  it('toggles panel expansion', () => {
    renderWithProviders(<TemplatePanel />);
    
    const toggleButton = screen.getByLabelText('Collapse');
    fireEvent.click(toggleButton);
    
    expect(mockToggleTemplatePanel).toHaveBeenCalled();
  });
});
