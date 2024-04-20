function fig = plotTriad(fig,o,M,scale,colorString)
    fig();
    quiver3(o(1),o(2),o(3), ...
        M(1,1),M(2,1),M(3,1), ...
        scale, ...
        'LineWidth',1, ...
        'Color',colorString)
    quiver3(o(1),o(2),o(3), ...
        M(1,2),M(2,2),M(3,2), ...
        scale, ...
        'LineWidth',1, ...
        'Color',colorString)
    quiver3(o(1),o(2),o(3), ...
        M(1,3),M(2,3),M(3,3), ...
        scale, ...
        'LineWidth',1, ...
        'Color',colorString)
end