library(data.table)
library(dplyr)
library(plyr)

Nrev = 25000

Direct = data.table(
  dgamma = rep(rep(seq(0,3e-4,length.out = 5), each = 5),5),
  x = rep(seq(-15e-3,15e-3,length.out = 5), times=5),
  Sx = rnorm(125, -.15, .1),
  Sy = rnorm(125, -.15, .1),
  Sz = rnorm(125, -.15, .1)
)

Reverse = data.table(
  dgamma = rep(rep(seq(0,3e-4,length.out = 5), each = 5),5),
  x = rep(seq(-15e-3,15e-3,length.out = 5), times=5),
  Sx = rnorm(125, -.15, .1),
  Sy = rnorm(125, -.15, .1),
  Sz = rnorm(125, -.15, .1)
)

Direct[,':='(
    ThXY = atan(Sy/Sx),
    ThXZ = atan(Sz/Sx),
    ThYZ = atan(Sy/Sz)
  )
]
Reverse[,':='(
    ThXY = atan(Sy/Sx),
    ThXZ = atan(Sz/Sx),
    ThYZ = atan(Sy/Sz)
  )
]
Direct[,':='(
  Wz = ThXY/Nrev,
  Wy = ThXZ/Nrev,
  Wx = ThYZ/Nrev
)]
Reverse[,':='(
  Wz = ThXY/Nrev,
  Wy = ThXZ/Nrev,
  Wx = ThYZ/Nrev
)]


names(Direct) <- paste0(names(Direct),".CW")
names(Reverse) <- paste0(names(Reverse),".CCW")

Direct[rep(1:nrow(Direct),each=nrow(Reverse)),]->Direct
Reverse[rep(1:nrow(Reverse),nrow(Reverse)),]->Reverse
cbind(Direct, Reverse) -> DR
DR[,':='(
    dWx = Wx.CW-Wx.CCW,
    dWy = Wy.CW-Wy.CCW,
    dWz = Wz.CW-Wz.CCW
)]
