import React from 'react'
import { Line } from 'react-chartjs-2'
import {
  Chart as ChartJS,
  CategoryScale,
  LinearScale,
  PointElement,
  LineElement,
  Title,
  Tooltip,
  Legend,
  Filler
} from 'chart.js'

ChartJS.register(
  CategoryScale,
  LinearScale,
  PointElement,
  LineElement,
  Title,
  Tooltip,
  Legend,
  Filler
)

interface IRPeak {
  wavenumber: number
  intensity: string
  group: string
  atoms?: string[]
}

interface IRGraphProps {
  peaks: IRPeak[]
  highlightPeak?: number
  onPeakHover?: (peak: IRPeak | null) => void
}

export default function IRGraph({ peaks, highlightPeak, onPeakHover }: IRGraphProps) {
  // Convert peaks to chart data
  const wavenumbers = peaks.map(p => p.wavenumber)
  const intensities = peaks.map(p => {
    const intensityMap: Record<string, number> = {
      'strong': 1.0,
      'medium': 0.6,
      'weak': 0.3,
      'broad': 0.8
    }
    return intensityMap[p.intensity.toLowerCase()] || 0.5
  })

  const data = {
    labels: wavenumbers,
    datasets: [
      {
        label: 'Transmittance',
        data: intensities.map(i => 100 - i * 100), // Invert for IR
        borderColor: highlightPeak !== undefined ? '#4676ff' : '#9aa0a6',
        backgroundColor: 'rgba(70, 118, 255, 0.1)',
        fill: true,
        tension: 0.4,
        pointRadius: 4,
        pointHoverRadius: 6,
      }
    ]
  }

  const options = {
    responsive: true,
    maintainAspectRatio: false,
    scales: {
      x: {
        reverse: true, // IR convention: high wavenumber on left
        title: {
          display: true,
          text: 'Wavenumber (cm⁻¹)'
        }
      },
      y: {
        min: 0,
        max: 100,
        title: {
          display: true,
          text: '% Transmittance'
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
            return `${peak.wavenumber} cm⁻¹ - ${peak.group} (${peak.intensity})`
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
      <Line data={data} options={options} />
    </div>
  )
}

