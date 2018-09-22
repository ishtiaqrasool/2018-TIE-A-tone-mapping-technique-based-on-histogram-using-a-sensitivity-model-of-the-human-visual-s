hdr = double(hdrread('memorial.hdr'));
ldr = ATT_TMO(hdr);
figure, imshow(ldr);
title('Adaptive TVI TMO by Khan et al.');