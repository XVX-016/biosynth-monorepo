import React from "react";
import { getTemplates, loadTemplate, placeTemplate } from "../../kernel/templateLoader";
import Button from "../ui/Button";

export const TemplatePicker: React.FC = () => {
  const templates = getTemplates();

  const handleLoad = (id: string) => {
    const tmpl = loadTemplate(id);
    if (!tmpl) return;

    placeTemplate(tmpl, { x: 0, y: 0, z: 0 });
  };

  return (
    <div className="p-4 bg-frostedGlass rounded-lg border border-chrome/30 w-full">
      <h3 className="text-lg font-bold mb-3 text-ivory">Templates</h3>

      <div className="grid grid-cols-2 gap-2">
        {templates.map((t) => (
          <Button key={t.id} variant="secondary" onClick={() => handleLoad(t.id)}>
            {t.name}
          </Button>
        ))}
      </div>
    </div>
  );
};

