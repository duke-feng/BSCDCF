function [res,w] = resp_new(h_1, h_2,kx)


responsef1 = fft2(response1);%%响应得傅里叶域
responsef2 = fft2(response2);
if numel(response2) ~= 1 && numel(response1) ~= 1
    sz = size(kx,1);
    rg = circshift(-floor((sz-1)/2):ceil((sz-1)/2), [0 -floor((sz-1)/2)]);
    cg = circshift(-floor((sz-1)/2):ceil((sz-1)/2), [0 -floor((sz-1)/2)]);
    [rs, cs] = ndgrid(rg,cg);
    d=1-exp(-1/sz*(rs.^2 + cs.^2));
    weight = 0.3:0.01:0.9;  %% 选取0.3到0.9 0.01为步长    
    L = zeros(numel(weight),1);
    i = 4;
    r1 = response1(:,:,i);
    r2 = response2(:,:,i);
    for j = 1:numel(weight)
        r = weight(j)*r1+(1-weight(j))*r2;%%手工特征与深度特征响应结合
        [r_max,id] = max(r(:));%%最大响应
        [max_id1,max_id2] = ind2sub(size(r),id);%%向量转换为矩阵
        r_r = (r_max - r)./circshift(r_max*d,[max_id1-1,max_id2-1]);%%循环移位
        ep = min(r_r(:));%%最小值
        L(j) = 0.1*(weight(j)^2+(1-weight(j))^2)-ep;
    end
    [~, jidx]=min(L(:));
    response = weight(jidx)*response1+(1-weight(jidx))*response2;
    responsef = weight(jidx)*responsef1+(1-weight(jidx))*responsef2; %%傅里叶域
    w = weight(jidx);
elseif numel(response2) ~= 1
    response = response2;
    responsef = responsef2;
    w = 0;
else
    response = response1;
    responsef = responsef1;
    w = 1;
end
end