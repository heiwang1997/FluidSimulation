pic_id = 0;
base_folder = '../1/';
for pic_id = 139
    rhoField = load([base_folder 'rho' int2str(pic_id)]);
    vxField = load([base_folder 'vx' int2str(pic_id)]);
    %vxWatch = vxField;
    vyField = load([base_folder 'vy' int2str(pic_id)]);
    % Interpolation.
    for i = 1:(size(vxField, 2)-1)
        vxField(:,i) = (vxField(:,i) + vxField(:,i+1)) ./ 2;
    end
    vxField = vxField(1:end, 1:end-1);
    for i = 1:(size(vyField, 1)-1)
        vyField(i,:) = (vyField(i,:) + vyField(i+1,:)) ./ 2;
    end
    vyField = vyField(1:end-1, 1:end);

    figure;
    quiver(vxField, vyField);
    title('速度场');
    saveas(gcf,['saved_fix_interface/vel' int2str(pic_id)],'bmp');

    figure;
    pcolor(rhoField);
    %quiver(vxWatch, zeros(size(vxWatch)));
    %quiver(vxField, vyField);
    title('密度场');
    colorbar;
    saveas(gcf,['saved_fix_interface/pic' int2str(pic_id)],'bmp');

    % close all;
end
