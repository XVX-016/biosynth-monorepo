import React, { useState, useRef, useEffect } from 'react';
import type { StudioMessage, StudioMode } from '../../types/studio';
import { StudioGateway } from '../../lib/ai/studioGateway';
import ChatMessage from './ChatMessage';

interface ChatInterfaceProps {
    mentorId: string;
    mode: StudioMode;
    onMessageSent?: (msg: StudioMessage) => void;
    onResponseReceived?: (msg: StudioMessage) => void;
}

export default function ChatInterface({ mentorId, mode, onMessageSent, onResponseReceived }: ChatInterfaceProps) {
    const [messages, setMessages] = useState<StudioMessage[]>([]);
    const [inputValue, setInputValue] = useState('');
    const [isLoading, setIsLoading] = useState(false);
    const scrollRef = useRef<HTMLDivElement>(null);

    useEffect(() => {
        if (scrollRef.current) {
            scrollRef.current.scrollTop = scrollRef.current.scrollHeight;
        }
    }, [messages]);

    const handleSubmit = async (e?: React.FormEvent) => {
        e?.preventDefault();
        if (!inputValue.trim() || isLoading) return;

        const userMessage: StudioMessage = {
            id: Date.now().toString(),
            role: 'user',
            content: inputValue,
            timestamp: new Date(),
        };

        setMessages(prev => [...prev, userMessage]);
        setInputValue('');
        setIsLoading(true);
        onMessageSent?.(userMessage);

        try {
            const response = await StudioGateway.chat({
                prompt: userMessage.content,
                mentorId,
                mode,
                history: messages,
            });

            let content = response.response;
            if (response.blocked) {
                content = "I'm sorry, but I can only assist with chemistry-related questions and molecular design. Is there a specific chemical structure or property you'd like to discuss?";
            }

            const mentorMessage: StudioMessage = {
                id: (Date.now() + 1).toString(),
                role: 'mentor',
                content: content,
                timestamp: new Date(),
                metadata: {
                    mode,
                }
            };

            setMessages(prev => [...prev, mentorMessage]);
            onResponseReceived?.(mentorMessage);
        } catch (error) {
            console.error("ChatInterface Error:", error);
        } finally {
            setIsLoading(false);
        }
    };

    return (
        <div className="flex flex-col h-full">
            {/* Chat Messages */}
            <div
                ref={scrollRef}
                className="flex-1 overflow-y-auto px-2 py-4 scrollbar-hide space-y-2"
            >
                {messages.length === 0 && (
                    <div className="flex flex-col items-center justify-center h-full opacity-40 text-center px-8">
                        <div className="w-12 h-12 border border-black/10 rounded-full flex items-center justify-center mb-4">
                            <svg width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
                                <path d="M21 15a2 2 0 0 1-2 2H7l-4 4V5a2 2 0 0 1 2-2h14a2 2 0 0 1 2 2z"></path>
                            </svg>
                        </div>
                        <p className="text-sm font-medium">No messages yet</p>
                        <p className="text-xs mt-1">Select a mentor and start designing.</p>
                    </div>
                )}
                {messages.map((msg) => (
                    <ChatMessage key={msg.id} message={msg} />
                ))}
                {isLoading && (
                    <div className="flex justify-start mb-4">
                        <div className="bg-transparent px-4 py-2">
                            <div className="flex gap-1">
                                <div className="w-1.5 h-1.5 bg-black/20 rounded-full animate-bounce" style={{ animationDelay: '0ms' }} />
                                <div className="w-1.5 h-1.5 bg-black/20 rounded-full animate-bounce" style={{ animationDelay: '150ms' }} />
                                <div className="w-1.5 h-1.5 bg-black/20 rounded-full animate-bounce" style={{ animationDelay: '300ms' }} />
                            </div>
                        </div>
                    </div>
                )}
            </div>

            {/* Input Area */}
            <form onSubmit={handleSubmit} className="relative mt-4">
                <input
                    type="text"
                    value={inputValue}
                    onChange={(e) => setInputValue(e.target.value)}
                    placeholder="Describe a molecule, reaction, or property..."
                    disabled={isLoading}
                    className="w-full bg-offwhite border border-lightGrey rounded-xl px-4 py-4 pr-12 text-black placeholder:text-midGrey focus:outline-none focus:ring-2 focus:ring-black/5 disabled:opacity-50"
                />
                <button
                    type="submit"
                    disabled={isLoading || !inputValue.trim()}
                    className="absolute right-3 top-1/2 -translate-y-1/2 bg-black text-white p-2 rounded-lg hover:bg-black/80 transition-colors disabled:bg-lightGrey disabled:cursor-not-allowed"
                >
                    <svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2.5" strokeLinecap="round" strokeLinejoin="round">
                        <line x1="22" y1="2" x2="11" y2="13"></line>
                        <polygon points="22 2 15 22 11 13 2 9 22 2"></polygon>
                    </svg>
                </button>
            </form>
        </div>
    );
}
