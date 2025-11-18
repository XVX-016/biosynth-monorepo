"""
Background job worker using RQ (Redis Queue)
"""
import os
import redis
from rq import Worker, Queue, Connection
from typing import Dict, Any
import json

# Redis connection
redis_url = os.getenv('REDIS_URL', 'redis://localhost:6379/0')
redis_conn = redis.from_url(redis_url)

# Queue names
QUEUE_NAMES = {
    'generation': 'molecule_generation',
    'prediction': 'property_prediction',
    'optimization': 'property_optimization',
    'reaction': 'reaction_processing',
}

# Create queues
queues = {
    name: Queue(name, connection=redis_conn)
    for name in QUEUE_NAMES.values()
}


def enqueue(task_name: str, payload: Dict[str, Any], queue_name: str = 'molecule_generation') -> str:
    """
    Enqueue a task
    
    Args:
        task_name: Name of the task function
        payload: Task payload
        queue_name: Queue to use
    
    Returns:
        Job ID
    """
    if queue_name not in QUEUE_NAMES.values():
        raise ValueError(f"Invalid queue name: {queue_name}")
    
    queue = queues[queue_name]
    job = queue.enqueue(
        task_name,
        **payload,
        job_timeout='10m',  # 10 minute timeout
        result_ttl=3600,  # Keep results for 1 hour
    )
    
    return job.id


def get_job_status(job_id: str) -> Dict[str, Any]:
    """
    Get job status
    
    Args:
        job_id: Job ID
    
    Returns:
        Job status dictionary
    """
    from rq.job import Job
    
    try:
        job = Job.fetch(job_id, connection=redis_conn)
        return {
            'id': job.id,
            'status': job.get_status(),
            'result': job.result,
            'error': str(job.exc_info) if job.exc_info else None,
            'created_at': job.created_at.isoformat() if job.created_at else None,
            'started_at': job.started_at.isoformat() if job.started_at else None,
            'ended_at': job.ended_at.isoformat() if job.ended_at else None,
        }
    except Exception as e:
        return {
            'id': job_id,
            'status': 'not_found',
            'error': str(e)
        }


# Task functions (these would be imported from actual modules)
def task_generate_molecule(prompt: str, constraints: Dict[str, Any] = None) -> Dict[str, Any]:
    """
    Generate molecule task (placeholder - implement with actual model)
    """
    # This would call the actual generation model
    # For now, return placeholder
    return {
        'smiles': 'C',
        'properties': {
            'stability': 0.5,
            'toxicity': 0.5,
            'solubility': 0.5,
            'bioavailability': 0.5,
            'novelty': 0.5,
        }
    }


def task_predict_properties(smiles: str) -> Dict[str, Any]:
    """
    Predict properties task
    """
    from backend.ai.predictor import get_predictor
    from backend.ai.featurizer import featurize_smiles
    
    try:
        features = featurize_smiles(smiles)
        if features is None:
            return {'error': 'Failed to featurize SMILES'}
        
        predictor = get_predictor()
        properties = predictor.predict(features)
        
        return {
            'smiles': smiles,
            'properties': properties
        }
    except Exception as e:
        return {'error': str(e)}


def task_optimize_properties(
    smiles: str,
    target_properties: Dict[str, float],
    iterations: int = 10
) -> Dict[str, Any]:
    """
    Optimize properties task (placeholder)
    """
    # This would implement property optimization
    return {
        'original_smiles': smiles,
        'optimized_smiles': smiles,
        'properties': target_properties,
        'iterations': iterations
    }


def task_process_reaction(
    reaction_type: str,
    reactants: str,
    **kwargs
) -> Dict[str, Any]:
    """
    Process reaction task
    """
    from backend.simulation.reaction_engine import reaction_engine
    
    try:
        if reaction_type == 'bond_break':
            result = reaction_engine.bond_break(
                reactants,
                kwargs.get('atom1_idx'),
                kwargs.get('atom2_idx')
            )
        elif reaction_type == 'bond_form':
            result = reaction_engine.bond_form(
                kwargs.get('smiles1', reactants),
                kwargs.get('smiles2'),
                kwargs.get('atom1_idx'),
                kwargs.get('atom2_idx')
            )
        elif reaction_type == 'ionize':
            result = reaction_engine.ionize(
                reactants,
                kwargs.get('atom_idx'),
                kwargs.get('charge', 1)
            )
        elif reaction_type == 'recombine':
            fragments = kwargs.get('fragments', [reactants])
            result = reaction_engine.recombine(fragments)
        else:
            return {'error': f'Unknown reaction type: {reaction_type}'}
        
        return {
            'success': result.success,
            'products': [
                {
                    'smiles': p.smiles,
                    'metadata': p.metadata
                }
                for p in result.products
            ],
            'error': result.error,
            'metadata': result.metadata
        }
    except Exception as e:
        return {'error': str(e)}


# Register task functions
TASK_FUNCTIONS = {
    'generate_molecule': task_generate_molecule,
    'predict_properties': task_predict_properties,
    'optimize_properties': task_optimize_properties,
    'process_reaction': task_process_reaction,
}


def start_worker(queue_names: list = None):
    """
    Start RQ worker
    
    Args:
        queue_names: List of queue names to listen to (None for all)
    """
    if queue_names is None:
        queue_names = list(QUEUE_NAMES.values())
    
    queues_to_listen = [queues[name] for name in queue_names if name in queues]
    
    with Connection(redis_conn):
        worker = Worker(queues_to_listen)
        worker.work()


if __name__ == '__main__':
    # Start worker
    import sys
    
    queue_names = sys.argv[1:] if len(sys.argv) > 1 else None
    start_worker(queue_names)

