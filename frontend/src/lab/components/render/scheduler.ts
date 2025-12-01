/**
 * Render Scheduler - RequestAnimationFrame queue and debounced tasks
 * 
 * Optimizes rendering performance by batching updates and throttling expensive operations.
 */

type Task = () => void;
type ScheduledTask = {
  task: Task;
  priority: 'high' | 'normal' | 'low';
  id: string;
};

export class RenderScheduler {
  private highPriorityQueue: ScheduledTask[] = [];
  private normalPriorityQueue: ScheduledTask[] = [];
  private lowPriorityQueue: ScheduledTask[] = [];
  private rafId: number | null = null;
  private isRunning = false;
  private debouncedTasks = new Map<string, { timeout: number; task: Task }>();

  /**
   * Schedule a high-priority task (rendering, user interactions)
   */
  scheduleHigh(task: Task, id?: string): string {
    const taskId = id || crypto.randomUUID();
    this.highPriorityQueue.push({ task, priority: 'high', id: taskId });
    this.start();
    return taskId;
  }

  /**
   * Schedule a normal-priority task
   */
  scheduleNormal(task: Task, id?: string): string {
    const taskId = id || crypto.randomUUID();
    this.normalPriorityQueue.push({ task, priority: 'normal', id: taskId });
    this.start();
    return taskId;
  }

  /**
   * Schedule a low-priority task (background work)
   */
  scheduleLow(task: Task, id?: string): string {
    const taskId = id || crypto.randomUUID();
    this.lowPriorityQueue.push({ task, priority: 'low', id: taskId });
    this.start();
    return taskId;
  }

  /**
   * Debounce a task (only execute after delay with no new calls)
   */
  debounce(id: string, task: Task, delay: number = 100): void {
    const existing = this.debouncedTasks.get(id);
    if (existing) {
      clearTimeout(existing.timeout);
    }

    const timeout = window.setTimeout(() => {
      task();
      this.debouncedTasks.delete(id);
    }, delay);

    this.debouncedTasks.set(id, { timeout, task });
  }

  /**
   * Throttle a task (execute at most once per interval)
   */
  throttle(id: string, task: Task, interval: number = 16): void {
    const existing = this.debouncedTasks.get(id);
    if (existing) {
      return; // Already scheduled
    }

    const timeout = window.setTimeout(() => {
      task();
      this.debouncedTasks.delete(id);
    }, interval);

    this.debouncedTasks.set(id, { timeout, task });
  }

  /**
   * Cancel a scheduled task
   */
  cancel(taskId: string): boolean {
    // Remove from queues
    const removeFromQueue = (queue: ScheduledTask[]) => {
      const index = queue.findIndex(t => t.id === taskId);
      if (index > -1) {
        queue.splice(index, 1);
        return true;
      }
      return false;
    };

    if (removeFromQueue(this.highPriorityQueue)) return true;
    if (removeFromQueue(this.normalPriorityQueue)) return true;
    if (removeFromQueue(this.lowPriorityQueue)) return true;

    // Cancel debounced task
    const debounced = this.debouncedTasks.get(taskId);
    if (debounced) {
      clearTimeout(debounced.timeout);
      this.debouncedTasks.delete(taskId);
      return true;
    }

    return false;
  }

  /**
   * Start the scheduler loop
   */
  private start(): void {
    if (this.isRunning) return;

    this.isRunning = true;
    this.rafId = requestAnimationFrame(() => this.tick());
  }

  /**
   * Scheduler tick
   */
  private tick(): void {
    // Execute high-priority tasks first
    while (this.highPriorityQueue.length > 0) {
      const task = this.highPriorityQueue.shift();
      if (task) {
        try {
          task.task();
        } catch (error) {
          console.error('Error executing high-priority task:', error);
        }
      }
    }

    // Execute normal-priority tasks (limit per frame)
    let normalCount = 0;
    const maxNormalPerFrame = 5;
    while (this.normalPriorityQueue.length > 0 && normalCount < maxNormalPerFrame) {
      const task = this.normalPriorityQueue.shift();
      if (task) {
        try {
          task.task();
        } catch (error) {
          console.error('Error executing normal-priority task:', error);
        }
        normalCount++;
      }
    }

    // Execute low-priority tasks (only if queues are empty)
    if (this.highPriorityQueue.length === 0 && this.normalPriorityQueue.length === 0) {
      if (this.lowPriorityQueue.length > 0) {
        const task = this.lowPriorityQueue.shift();
        if (task) {
          try {
            task.task();
          } catch (error) {
            console.error('Error executing low-priority task:', error);
          }
        }
      }
    }

    // Continue if there are more tasks
    if (this.hasPendingTasks()) {
      this.rafId = requestAnimationFrame(() => this.tick());
    } else {
      this.isRunning = false;
      this.rafId = null;
    }
  }

  /**
   * Check if there are pending tasks
   */
  private hasPendingTasks(): boolean {
    return (
      this.highPriorityQueue.length > 0 ||
      this.normalPriorityQueue.length > 0 ||
      this.lowPriorityQueue.length > 0
    );
  }

  /**
   * Clear all queues
   */
  clear(): void {
    this.highPriorityQueue = [];
    this.normalPriorityQueue = [];
    this.lowPriorityQueue = [];

    // Clear debounced tasks
    for (const { timeout } of this.debouncedTasks.values()) {
      clearTimeout(timeout);
    }
    this.debouncedTasks.clear();

    // Cancel RAF
    if (this.rafId !== null) {
      cancelAnimationFrame(this.rafId);
      this.rafId = null;
    }

    this.isRunning = false;
  }

  /**
   * Get queue statistics
   */
  getStats(): {
    high: number;
    normal: number;
    low: number;
    debounced: number;
  } {
    return {
      high: this.highPriorityQueue.length,
      normal: this.normalPriorityQueue.length,
      low: this.lowPriorityQueue.length,
      debounced: this.debouncedTasks.size,
    };
  }
}

// Singleton instance
export const renderScheduler = new RenderScheduler();

