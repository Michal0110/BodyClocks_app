// Custom tooltip handling for visNetwork
$(document).ready(function() {
  const EDGE_PADDING = 10;
  const POINTER_OFFSET = 12;
  const FLIP_START_FRACTION = 2 / 3;
  const MAX_TOOLTIP_WIDTH = 400;
  const MIN_TOOLTIP_WIDTH = 180;

  let pointer = {
    x: window.innerWidth / 2,
    y: window.innerHeight / 2
  };
  let pendingPosition = false;

  function clamp(value, min, max) {
    if (max < min) {
      return min;
    }
    return Math.min(Math.max(value, min), max);
  }

  function setStyle(element, property, value) {
    if (element.style[property] !== value) {
      element.style[property] = value;
    }
  }

  function getGraphContainer(element) {
    return element.closest('[id^="graph"]');
  }

  function getGraphRect(element) {
    const graph = getGraphContainer(element);

    if (!graph) {
      return {
        left: 0,
        right: window.innerWidth,
        top: 0,
        bottom: window.innerHeight,
        width: window.innerWidth,
        height: window.innerHeight
      };
    }

    const rect = graph.getBoundingClientRect();

    return {
      left: rect.left,
      right: rect.right,
      top: rect.top,
      bottom: rect.bottom,
      width: rect.width,
      height: rect.height
    };
  }

  function getMaxTooltipWidth(element) {
    const graphRect = getGraphRect(element);
    const viewportWidth = Math.max(80, window.innerWidth - (EDGE_PADDING * 2));
    const graphWidth = Math.max(80, graphRect.width - (EDGE_PADDING * 2));
    const availableWidth = Math.min(
      viewportWidth,
      Math.max(MIN_TOOLTIP_WIDTH, graphWidth)
    );

    return Math.min(MAX_TOOLTIP_WIDTH, availableWidth);
  }

  function isVisNetworkWrapperTooltip(element) {
    const parent = element.parentElement;

    if (!parent || !parent.id || !parent.id.startsWith('graph')) {
      return false;
    }

    const inlineStyle = (element.getAttribute('style') || '').toLowerCase();
    const computedStyle = window.getComputedStyle(element);

    return computedStyle.position === 'fixed' &&
      inlineStyle.indexOf('visibility') !== -1 &&
      !element.classList.contains('vis-network');
  }

  function getTooltipElements() {
    return $('.vis-tooltip, [id^="graph"] > div').filter(function() {
      return this.classList.contains('vis-tooltip') ||
        isVisNetworkWrapperTooltip(this);
    });
  }

  function normaliseTooltipSize(tooltip) {
    const maxWidth = getMaxTooltipWidth(tooltip);

    setStyle(tooltip, 'boxSizing', 'border-box');
    setStyle(tooltip, 'right', 'auto');
    setStyle(tooltip, 'width', 'max-content');
    setStyle(tooltip, 'maxWidth', maxWidth + 'px');
    setStyle(tooltip, 'whiteSpace', 'normal');
    setStyle(tooltip, 'wordBreak', 'normal');
    setStyle(tooltip, 'wordWrap', 'break-word');
    setStyle(tooltip, 'overflowWrap', 'break-word');
  }

  function positionTooltip(tooltip) {
    const computedStyle = window.getComputedStyle(tooltip);

    if (computedStyle.visibility === 'hidden' || computedStyle.display === 'none') {
      normaliseTooltipSize(tooltip);
      return;
    }

    normaliseTooltipSize(tooltip);

    const graphRect = getGraphRect(tooltip);
    const tooltipRect = tooltip.getBoundingClientRect();
    const tooltipWidth = tooltipRect.width;
    const tooltipHeight = tooltipRect.height;
    const viewportLeft = EDGE_PADDING;
    const viewportRight = window.innerWidth - EDGE_PADDING;
    const viewportTop = EDGE_PADDING;
    const viewportBottom = window.innerHeight - EDGE_PADDING;
    const graphPointerX = pointer.x - graphRect.left;
    const preferLeft = graphPointerX >= graphRect.width * FLIP_START_FRACTION;
    const roomRight = viewportRight - pointer.x - POINTER_OFFSET;
    const roomLeft = pointer.x - viewportLeft - POINTER_OFFSET;
    const flipLeft = preferLeft || (roomRight < tooltipWidth && roomLeft > roomRight);
    let left = flipLeft
      ? pointer.x - tooltipWidth - POINTER_OFFSET
      : pointer.x + POINTER_OFFSET;
    let top = pointer.y - 20;

    left = clamp(left, viewportLeft, viewportRight - tooltipWidth);
    top = clamp(top, viewportTop, viewportBottom - tooltipHeight);

    setStyle(tooltip, 'left', left + 'px');
    setStyle(tooltip, 'top', top + 'px');
  }

  function positionVisibleTooltips() {
    getTooltipElements().each(function() {
      positionTooltip(this);
    });
  }

  function schedulePositionVisibleTooltips() {
    if (pendingPosition) {
      return;
    }

    pendingPosition = true;
    window.requestAnimationFrame(function() {
      pendingPosition = false;
      positionVisibleTooltips();
    });
  }

  $(document).on('mousemove', function(event) {
    pointer.x = event.clientX;
    pointer.y = event.clientY;
    schedulePositionVisibleTooltips();
  });

  $(document).on('mouseover', '.vis-tooltip, [id^="graph"] > div', function() {
    schedulePositionVisibleTooltips();
  });

  $(window).on('resize', schedulePositionVisibleTooltips);

  const tooltipObserver = new MutationObserver(schedulePositionVisibleTooltips);

  tooltipObserver.observe(document.body, {
    childList: true,
    subtree: true,
    attributes: true,
    attributeFilter: ['style', 'class']
  });
});
