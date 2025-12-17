import React from 'react';
import { MENTORS } from '../../config/mentors';
import StudioMentorCard from './StudioMentorCard';

interface MentorStripProps {
    activeMentorId?: string;
    onMentorSelect: (mentorId: string) => void;
}

export default function MentorStrip({ activeMentorId, onMentorSelect }: MentorStripProps) {
    return (
        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-6">
            {MENTORS.map((mentor) => (
                <StudioMentorCard
                    key={mentor.id}
                    mentor={mentor}
                    isActive={activeMentorId === mentor.id}
                    onClick={onMentorSelect}
                />
            ))}
        </div>
    );
}
