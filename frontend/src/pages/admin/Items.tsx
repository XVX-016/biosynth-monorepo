import React, { useState, useEffect } from 'react';
import { motion } from 'framer-motion';
import Card from '../../components/ui/Card';
import Button from '../../components/ui/Button';
import { listItems, createItem, updateItem, deleteItem } from '../../lib/api';
import type { Item } from '../../lib/api';

export default function AdminItems() {
  const [items, setItems] = useState<Item[]>([]);
  const [loading, setLoading] = useState(false);
  const [showForm, setShowForm] = useState(false);
  const [editingItem, setEditingItem] = useState<Item | null>(null);
  const [formData, setFormData] = useState({
    name: '',
    smiles: '',
    description: '',
    tags: [] as string[],
    status: 'in-stock' as 'in-stock' | 'sold-out' | 'archived',
    stock: 0,
    structure_file: null as File | null,
  });
  const [tagInput, setTagInput] = useState('');

  useEffect(() => {
    loadItems();
  }, []);

  const loadItems = async () => {
    setLoading(true);
    try {
      const data = await listItems();
      setItems(data);
    } catch (error) {
      console.error('Failed to load items:', error);
    } finally {
      setLoading(false);
    }
  };

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    setLoading(true);
    try {
      if (editingItem) {
        await updateItem(editingItem.id, formData);
      } else {
        await createItem(formData);
      }
      setShowForm(false);
      setEditingItem(null);
      resetForm();
      loadItems();
    } catch (error) {
      console.error('Failed to save item:', error);
      alert('Failed to save item');
    } finally {
      setLoading(false);
    }
  };

  const handleEdit = (item: Item) => {
    setEditingItem(item);
    setFormData({
      name: item.name,
      smiles: item.smiles || '',
      description: item.description || '',
      tags: item.tags || [],
      status: item.status,
      stock: item.stock || 0,
      structure_file: null,
    });
    setShowForm(true);
  };

  const handleDelete = async (id: number) => {
    if (!confirm('Are you sure you want to delete this item?')) return;
    
    setLoading(true);
    try {
      await deleteItem(id);
      loadItems();
    } catch (error) {
      console.error('Failed to delete item:', error);
      alert('Failed to delete item');
    } finally {
      setLoading(false);
    }
  };

  const resetForm = () => {
    setFormData({
      name: '',
      smiles: '',
      description: '',
      tags: [],
      status: 'in-stock',
      stock: 0,
      structure_file: null,
    });
    setTagInput('');
  };

  const addTag = () => {
    if (tagInput.trim() && !formData.tags.includes(tagInput.trim())) {
      setFormData({
        ...formData,
        tags: [...formData.tags, tagInput.trim()],
      });
      setTagInput('');
    }
  };

  const removeTag = (tag: string) => {
    setFormData({
      ...formData,
      tags: formData.tags.filter((t) => t !== tag),
    });
  };

  const handleFileChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (file) {
      setFormData({ ...formData, structure_file: file });
    }
  };

  return (
    <motion.div
      initial={{ opacity: 0 }}
      animate={{ opacity: 1 }}
      transition={{ duration: 0.3 }}
      className="space-y-6"
    >
      <div className="flex items-center justify-between">
        <h1 className="text-3xl font-bold text-ivory">Admin - Items</h1>
        <Button onClick={() => {
          resetForm();
          setEditingItem(null);
          setShowForm(true);
        }}>
          Add New Item
        </Button>
      </div>

      {showForm && (
        <Card>
          <form onSubmit={handleSubmit} className="space-y-4">
            <div>
              <label className="block text-sm font-medium mb-1 text-chrome">Name</label>
              <input
                type="text"
                value={formData.name}
                onChange={(e) => setFormData({ ...formData, name: e.target.value })}
                className="w-full px-3 py-2 border border-chrome/20 bg-frostedGlass text-ivory rounded-lg focus:ring-2 focus:ring-neonCyan/50 focus:border-neonCyan/50 outline-none"
                required
              />
            </div>

            <div>
              <label className="block text-sm font-medium mb-1 text-chrome">SMILES</label>
              <input
                type="text"
                value={formData.smiles}
                onChange={(e) => setFormData({ ...formData, smiles: e.target.value })}
                className="w-full px-3 py-2 border border-chrome/20 bg-frostedGlass text-ivory rounded-lg focus:ring-2 focus:ring-neonCyan/50 focus:border-neonCyan/50 outline-none"
              />
            </div>

            <div>
              <label className="block text-sm font-medium mb-1 text-chrome">Description</label>
              <textarea
                value={formData.description}
                onChange={(e) => setFormData({ ...formData, description: e.target.value })}
                className="w-full px-3 py-2 border border-chrome/20 bg-frostedGlass text-ivory rounded-lg focus:ring-2 focus:ring-neonCyan/50 focus:border-neonCyan/50 outline-none"
                rows={3}
              />
            </div>

            <div>
              <label className="block text-sm font-medium mb-1 text-chrome">Tags</label>
              <div className="flex gap-2 mb-2">
                <input
                  type="text"
                  value={tagInput}
                  onChange={(e) => setTagInput(e.target.value)}
                  onKeyPress={(e) => {
                    if (e.key === 'Enter') {
                      e.preventDefault();
                      addTag();
                    }
                  }}
                  className="flex-1 px-3 py-2 border border-chrome/20 bg-frostedGlass text-ivory rounded-lg focus:ring-2 focus:ring-neonCyan/50 focus:border-neonCyan/50 outline-none placeholder:text-chrome"
                  placeholder="Add tag..."
                />
                <Button type="button" onClick={addTag}>Add</Button>
              </div>
              <div className="flex flex-wrap gap-2">
                {formData.tags.map((tag) => (
                  <span
                    key={tag}
                    className="px-2 py-1 bg-frostedGlass border border-chrome/20 rounded text-sm flex items-center gap-1 text-chrome"
                  >
                    {tag}
                    <button
                      type="button"
                      onClick={() => removeTag(tag)}
                      className="text-chrome/70 hover:text-ivory transition-colors"
                    >
                      ×
                    </button>
                  </span>
                ))}
              </div>
            </div>

            <div className="grid grid-cols-2 gap-4">
              <div>
                <label className="block text-sm font-medium mb-1 text-chrome">Status</label>
                <select
                  value={formData.status}
                  onChange={(e) => setFormData({ ...formData, status: e.target.value as any })}
                  className="w-full px-3 py-2 border border-chrome/20 bg-frostedGlass text-ivory rounded-lg focus:ring-2 focus:ring-neonCyan/50 focus:border-neonCyan/50 outline-none"
                >
                  <option value="in-stock">In Stock</option>
                  <option value="sold-out">Sold Out</option>
                  <option value="archived">Archived</option>
                </select>
              </div>

              <div>
                <label className="block text-sm font-medium mb-1 text-chrome">Stock</label>
                <input
                  type="number"
                  value={formData.stock}
                  onChange={(e) => setFormData({ ...formData, stock: parseInt(e.target.value) || 0 })}
                  className="w-full px-3 py-2 border border-chrome/20 bg-frostedGlass text-ivory rounded-lg focus:ring-2 focus:ring-neonCyan/50 focus:border-neonCyan/50 outline-none"
                  min="0"
                />
              </div>
            </div>

            <div>
              <label className="block text-sm font-medium mb-1 text-chrome">Structure File</label>
              <input
                type="file"
                onChange={handleFileChange}
                accept=".mol,.sdf,.pdb"
                className="w-full px-3 py-2 border border-chrome/20 bg-frostedGlass text-ivory rounded-lg focus:ring-2 focus:ring-neonCyan/50 focus:border-neonCyan/50 outline-none file:mr-4 file:py-1 file:px-3 file:rounded file:border-0 file:text-sm file:bg-plasma-neon file:text-ionBlack file:cursor-pointer"
              />
            </div>

            <div className="flex gap-2">
              <Button type="submit" disabled={loading}>
                {editingItem ? 'Update' : 'Create'} Item
              </Button>
              <Button
                type="button"
                onClick={() => {
                  setShowForm(false);
                  setEditingItem(null);
                  resetForm();
                }}
              >
                Cancel
              </Button>
            </div>
          </form>
        </Card>
      )}

      <Card>
        <div className="space-y-4">
          {loading && items.length === 0 ? (
            <div className="text-center py-8 text-chrome">Loading...</div>
          ) : items.length === 0 ? (
            <div className="text-center py-8 text-chrome">No items found</div>
          ) : (
            <div className="overflow-x-auto">
              <table className="w-full">
                <thead>
                  <tr className="border-b border-chrome/20">
                    <th className="text-left py-2 px-4 text-chrome">Name</th>
                    <th className="text-left py-2 px-4 text-chrome">SMILES</th>
                    <th className="text-left py-2 px-4 text-chrome">Status</th>
                    <th className="text-left py-2 px-4 text-chrome">Stock</th>
                    <th className="text-left py-2 px-4 text-chrome">Tags</th>
                    <th className="text-right py-2 px-4 text-chrome">Actions</th>
                  </tr>
                </thead>
                <tbody>
                  {items.map((item) => (
                    <tr key={item.id} className="border-b border-chrome/10">
                      <td className="py-2 px-4 text-ivory">{item.name}</td>
                      <td className="py-2 px-4 font-mono text-sm text-chrome">{item.smiles || '—'}</td>
                      <td className="py-2 px-4">
                        <span
                          className={`px-2 py-1 rounded text-xs ${
                            item.status === 'in-stock'
                              ? 'bg-plasmaTeal/20 text-plasmaTeal border border-plasmaTeal/30'
                              : item.status === 'sold-out'
                              ? 'bg-violetEdge/20 text-violetEdge border border-violetEdge/30'
                              : 'bg-chrome/20 text-chrome border border-chrome/30'
                          }`}
                        >
                          {item.status}
                        </span>
                      </td>
                      <td className="py-2 px-4 text-ivory">{item.stock || 0}</td>
                      <td className="py-2 px-4">
                        <div className="flex flex-wrap gap-1">
                          {item.tags?.slice(0, 3).map((tag) => (
                            <span
                              key={tag}
                              className="px-1 py-0.5 bg-frostedGlass border border-chrome/20 rounded text-xs text-chrome"
                            >
                              {tag}
                            </span>
                          ))}
                          {item.tags && item.tags.length > 3 && (
                            <span className="text-xs text-chrome/70">+{item.tags.length - 3}</span>
                          )}
                        </div>
                      </td>
                      <td className="py-2 px-4 text-right">
                        <div className="flex gap-2 justify-end">
                          <Button
                            size="sm"
                            onClick={() => handleEdit(item)}
                          >
                            Edit
                          </Button>
                          <Button
                            size="sm"
                            variant="danger"
                            onClick={() => handleDelete(item.id)}
                          >
                            Delete
                          </Button>
                        </div>
                      </td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          )}
        </div>
      </Card>
    </motion.div>
  );
}

