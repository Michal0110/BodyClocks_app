$(document).ready(function() {
    // Create SVG element
    const container = document.getElementById('sine-wave-container');
    const svgns = "http://www.w3.org/2000/svg";
    const svg = document.createElementNS(svgns, "svg");
    
    // Configuration
    const waves = [
        { frequency: 4, speed: 0.02, amplitude: 0.5, color: '#bdc3c7', opacity: 0.4, width: 3 },
        { frequency: 3.5, speed: -0.015, amplitude: 0.45, color: '#bdc3c7', opacity: 0.3, width: 3 },
        { frequency: 4.5, speed: 0.025, amplitude: 0.4, color: '#bdc3c7', opacity: 0.35, width: 3 }
    ];
    
    // Animation variables
    let animationFrameId;
    let phases = [0, Math.PI * 0.5, Math.PI];
    
    // Function to generate a single sine wave path
    function generateSineWave(width, height, waveConfig, phase) {
        const amplitude = height * waveConfig.amplitude;
        const frequency = waveConfig.frequency;
        const points = [];
        
        // Use more points for smoother curves
        const steps = Math.max(100, width);
        
        for (let i = 0; i <= steps; i++) {
            const x = (i / steps) * width;
            const wavePhase = (x / width) * Math.PI * frequency + phase;
            const y = amplitude * Math.sin(wavePhase) + (height / 2);
            points.push(`${x},${y}`);
        }
        
        return points.join(' ');
    }
    
    // Create path elements for each wave
    const paths = waves.map((_, index) => {
        const path = document.createElementNS(svgns, "polyline");
        path.setAttribute('fill', 'none');
        path.setAttribute('stroke', waves[index].color);
        path.setAttribute('stroke-width', waves[index].width);
        path.setAttribute('stroke-opacity', waves[index].opacity);
        svg.appendChild(path);
        return path;
    });
    
    // Animation function
    function animate() {
        // Update phases
        phases = phases.map((phase, index) => {
            return (phase + waves[index].speed) % (Math.PI * 2);
        });
        
        // Update paths
        const width = container.offsetWidth;
        const height = container.offsetHeight;
        
        paths.forEach((path, index) => {
            const points = generateSineWave(width, height, waves[index], phases[index]);
            path.setAttribute('points', points);
        });
        
        // Continue animation loop
        animationFrameId = requestAnimationFrame(animate);
    }
    
    // Function to set up SVG and start animation
    function setupAnimation() {
        const width = container.offsetWidth;
        const height = container.offsetHeight;
        
        // Set SVG attributes
        svg.setAttribute('width', width);
        svg.setAttribute('height', height);
        svg.setAttribute('viewBox', `0 0 ${width} ${height}`);
        
        // Ensure container has the SVG
        if (!container.contains(svg)) {
            container.appendChild(svg);
        }
        
        // Start animation loop
        if (animationFrameId) {
            cancelAnimationFrame(animationFrameId);
        }
        animationFrameId = requestAnimationFrame(animate);
    }
    
    // Initial setup
    setupAnimation();
    
    // Update on window resize
    let resizeTimeout;
    window.addEventListener('resize', function() {
        clearTimeout(resizeTimeout);
        resizeTimeout = setTimeout(setupAnimation, 100);
    });
    
    // Clean up on page leave
    window.addEventListener('beforeunload', function() {
        if (animationFrameId) {
            cancelAnimationFrame(animationFrameId);
        }
    });
});