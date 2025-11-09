import { describe, it, expect } from 'vitest'
import { render } from '@testing-library/react'
import MoleculeViewer from '../MoleculeViewer'

describe('MoleculeViewer', () => {
  it('renders without crashing', () => {
    const { container } = render(<MoleculeViewer />)
    expect(container).toBeTruthy()
  })
})

