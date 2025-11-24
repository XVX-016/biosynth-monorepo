import React from 'react'
import { Bar } from 'react-chartjs-2'
import {
  Chart as ChartJS,
  CategoryScale,
  LinearScale,
  BarElement,
  Title,
  Tooltip,
  Legend
} from 'chart.js'

ChartJS.register(
  CategoryScale,
  LinearScale,
  BarElement,
  Title,
  Tooltip,
  Legend
)

interface MassPeak {
  'm/z': number
  intensity: number
  fragment_type?: string
  fragment_atoms?: string[]
}

interface MassGraphProps {
  peaks: MassPeak[]
  molecularIon: number
  highlightPeak?: number
  onPeakHover?: (peak: MassPeak | null) => void
}

export default function MassGraph({ peaks, molecularIon, highlightPeak, onPeakHover }: MassGraphProps) {
  const mzValues = peaks.map(p => p['m/z'])
  const intensities = peaks.map(p => p.intensity * 100) // Scale to percentage

  const data = {
    labels: mzValues.map(m => `${m.toFixed(1)}`),
    datasets: [
      {
        label: 'Relative Intensity (%)',
        data: intensities,
        backgroundColor: highlightPeak !== undefined
          ? mzValues.map((_, i) => i === highlightPeak ? '#4676ff' : '#9aa0a6')
          : '#9aa0a6',
        borderColor: '#4676ff',
        borderWidth: 1
      }
    ]
  }

  const options = {
    responsive: true,
    maintainAspectRatio: false,
    scales: {
      x: {
        title: {
          display: true,
          text: 'm/z'
        }
      },
      y: {
        min: 0,
        max: 100,
        title: {
          display: true,
          text: 'Relative Intensity (%)'
        }
      }
    },
    plugins: {
      legend: {
        display: false
      },
      tooltip: {
        callbacks: {
          label: (context: any) => {
            const peak = peaks[context.dataIndex]
            let label = `m/z ${peak['m/z'].toFixed(1)} - ${(peak.intensity * 100).toFixed(1)}%`
            if (peak.fragment_type) {
              label += ` (${peak.fragment_type})`
            }
            return label
          }
        }
      }
    },
    onHover: (event: any, elements: any[]) => {
      if (elements.length > 0) {
        const index = elements[0].index
        onPeakHover?.(peaks[index])
      } else {
        onPeakHover?.(null)
      }
    }
  }

  return (
    <div className="w-full h-64">
      <Bar data={data} options={options} />
      <div className="mt-2 text-xs text-gray-600 text-center">
        Molecular Ion: {molecularIon.toFixed(1)} m/z
      </div>
    </div>
  )
}

