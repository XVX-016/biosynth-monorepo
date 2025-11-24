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

interface NMRPeak {
  shift: number
  multiplicity?: string
  integration?: number
  atom?: string
}

interface NMRGraphProps {
  peaks: NMRPeak[]
  type: '1H' | '13C'
  highlightPeak?: number
  onPeakHover?: (peak: NMRPeak | null) => void
}

export default function NMRGraph({ peaks, type, highlightPeak, onPeakHover }: NMRGraphProps) {
  const shifts = peaks.map(p => p.shift)
  const integrations = peaks.map(p => p.integration || 1)

  const data = {
    labels: shifts.map(s => `${s} ppm`),
    datasets: [
      {
        label: type === '1H' ? 'Integration' : 'Intensity',
        data: integrations,
        backgroundColor: highlightPeak !== undefined 
          ? shifts.map((_, i) => i === highlightPeak ? '#4676ff' : '#9aa0a6')
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
        reverse: true, // NMR convention: high ppm on left
        title: {
          display: true,
          text: 'Chemical Shift (ppm)'
        }
      },
      y: {
        title: {
          display: true,
          text: type === '1H' ? 'Integration' : 'Intensity'
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
            let label = `${peak.shift} ppm`
            if (peak.multiplicity) {
              label += ` (${peak.multiplicity})`
            }
            if (peak.integration) {
              label += ` - ${peak.integration}H`
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
    </div>
  )
}

