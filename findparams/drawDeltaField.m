pic_id = 10;
rhoField = load(['../1/rho' int2str(pic_id)]);
vxField = load(['../1/vx' int2str(pic_id)]);
vyField = load(['../1/vy' int2str(pic_id)]);
vxField = vxField(1:end-1, 1:end);
vyField = vyField(1:end, 1:end-1);

quiver(vxField, vyField);
title('速度变化场（相对上一个dt）');

figure;
pcolor(rhoField);
title('密度变化场（相对上一个dt）');
colorbar;
