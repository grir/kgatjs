"use strict";
/*
Mokomoji js bilioteka KGAT kursui
Autorius: I.Grinis
*/

////////////////////////////////////////////////////
/// class MatNxN
////////////////////////////////////////////////////


class MatNxN { //row-major order !!! (Irus Grinis) 
  constructor(dim, diag_val) {
    this.array = new Array(dim);
    this.dim = dim;
    for(var i=0;i<dim;i++)
       for(var j=0;j<dim;j++){
          var addr = i * this.dim + j;  
          var v = 0.0; 
          if (i==j) v = diag_val;
          this.array[addr] = v;
       }    
  }

  getAt(i,j){
      var addr = i * this.dim + j;  
      return this.array[addr];
  }
  
  
  getSubMatrix(i,j){
      var subm = [];
      var n = this.dim;
      for(var a=0; a<n; a++)
         for(var b=0; b<n; b++)
            if ((i!==a)&&(j!==b))
               subm.push(this.getAt(a,b));
      var res = new MatNxN(n-1,0);         
      res.setByRows(subm);
      return res;
  }
  
  
  setByRows(array){
    var n = this.dim;
    n = n * n;
    for(var i=0;i<n;i++)
      this.array[i]=array[i];
  }
  
  setByColumns(array){
    var n = this.dim;
    var k = 0;
    for(var i=0;i<n;i++)
       for(var j=0;j<n;j++){
         this.setAt(j, i, array[k]);
         k = k + 1;
       } 
  }
    
  setAt(i,j, val){
      var addr = i * this.dim + j;  
      this.array[addr] = val;
    
  }
  
  getTranspose(){
    var n = this.dim;
    var tM = new MatNxN(n,0);
    for(var i=0;i<n;i++)
       for(var j=0;j<n;j++){
         tM.setAt(j, i, this.getAt(i,j));
       } 
   return tM;      
  }
  
  getArrayByRows(){
    return this.array.slice();   
  }
  
  getArrayByColumns(){
    return this.getTranspose().array.slice();   
  }
  
  
  toString(tofix=3){
    var s="";
    var n = this.dim;
    for(var i=0;i<n;i++){
      for(var j=0;j<n;j++){
        s = s + this.getAt(i,j).toFixed(tofix) + " ";  
      }
      s = s + "\n";
    }
    return s;        
  }
  
  
  static getDet(mA){
    var n = mA.dim;
    if(mA.dim <= 3){
      if(mA.dim === 3){
         var mA3 = new Mat3x3(0);
         mA3.setByRows(mA.array);
         return mA3.det();
      }
      if(mA.dim === 2){
         var mA2 = new Mat2x2(0);
         mA2.setByRows(mA.array);
         return mA2.det();
      }
      else return mA.array[0]; 
    }
    else{
      var ml = 1.0;
      var sum = 0.0;
      for(var j = 0; j < n; j++){
        var A0j = mA.getSubMatrix(0,j);
        var dt = MatNxN.getDet(A0j)
        //console.log("Call for \n" + A0j.toString() + "\n det <- " + dt);
        sum = sum + ml * dt * mA.getAt(0,j);
        ml = -ml;
      }
      return sum;
    }    
  }
  
  static getCof(mA, i, j){
    var subM = mA.getSubMatrix(i,j);
    var minor = MatNxN.getDet(subM);
    var sn = 1.0;
    if(((i+j)%2)!==0) sn = -1.0;
    return minor * sn;
  }
  
////////////////////////////////////////////////////
//////  mA^(-1)
////
//  

  static getInverse(mA, i, j){
    var inv =[];
    var n = mA.dim;
    var dt = MatNxN.getDet(mA);
    if (dt === 0) return null;
    for(var a=0;a<n;a++)
       for(var b=0;b<n;b++){
         var cf = MatNxN.getCof(mA, b, a);
         inv.push(cf/dt);
       }
          
    var invM = new MatNxN(n,0);
    invM.setByRows(inv);
    return invM;
  }
  

////////////////////////////////////////////////////
//////  sum <- mA + mB
////
//  
  static add(mA,mB){
      var sum = new MatNxN(mA.dim, 0);
      var n = mA.dim;
      var n2 = n * n;
      for(var i=0;i<n2;i++)
        sum.array[i]=mA.array[i]+mB.array[i];
      return sum; 
  }

////////////////////////////////////////////////////
//////  mil <- mA * mB
////
//  
  
  static mMul(mA,mB){
      var mul = new MatNxN(mA.dim, 0);
      var n = mA.dim;
      for(var i=0;i<n;i++)
        for(var j=0;j<n;j++){
           var s = 0.0;
           for(var k=0;k<n;k++)
              s = s + mA.getAt(i,k) * mB.getAt(k,j);
           mul.setAt(i,j, s);    
        }      
      return mul; 
    }
}

////////////////////////////////////////////////////
/// class Mat2x2
////////////////////////////////////////////////////

class Mat2x2 extends MatNxN{
   constructor(diag_val) {
     super(2, diag_val);
   }
   
   det(){
     var a = this.array;
     return  a[0]*a[3] - a[2]*a[1];
   }
   
}

////////////////////////////////////////////////////
/// class Mat3x3
////////////////////////////////////////////////////



class Mat3x3 extends MatNxN{
   constructor(diag_val) {
     super(3, diag_val);
   }
   det(){
     var a = this.array;
     return  a[0]*a[4]*a[8] + a[6]*a[1]*a[5] + a[2]*a[3]*a[7] -
             a[6]*a[4]*a[2] - a[1]*a[3]*a[8] - a[0]*a[7]*a[5];
   }
   cof(i,j){
     
     
   }

}

////////////////////////////////////////////////////
/// class Mat4x4
////////////////////////////////////////////////////


class Mat4x4 extends MatNxN{
   constructor(diag_val) {
     super(4, diag_val);
   }
   //Translation
   setTranslation(x,y,z){
     var a = this.array;
     a[0]=1.0;a[5]=1.0;a[10]=1.0;a[15]=1.0;
     a[3]=x;a[7]=y;a[11]=z;
   }
   //Scale
   setScale(x,y,z){
     var a = this.array;
     a[0]=x;a[5]=y;a[10]=z;a[15]=1.0;
   }
   //Rotation
   setRotation(angleDegrees,x,y,z){
     var v3 = new Vec3(x,y,z);
     var n = v3.norm();
     var I = new Mat4x4(1);
     
     if (n===0.0){ this.setByRows(I.array); return; }
     
     v3.normalize();
     var theta = angleDegrees * (Math.PI / 180);
     var cs = Math.cos(theta);
     var cs1 = 1 - cs;
     var sn = Math.sin(theta);
     ///////////////////////////////////////
     var Ucross = new Mat4x4(0);
     var ua = Ucross.array;
     var v = v3.array;
     ua[1] = -v[2]; ua[2] = -v[1];
     ua[4] = v[2];  ua[6] = -v[0];
     ua[8] = -v[1];  ua[9] = v[0];
     for(var i=0;i<12;i++)
        ua[i] = sn * ua[i];
     //ua[15] = 1.0;
   
     
     var Utensor = new Mat4x4(0);
     for(i=0;i<3;i++)
       for(var j=0;j<3;j++)
          Utensor.setAt(i,j,cs1 * v[i] * v[j]);
          

     I.setScale(cs,cs,cs); // cos(theta) * I 

     var result = Mat4x4.add(Mat4x4.add(I,Ucross),Utensor);
     this.setByRows(result.array);
   }

   setOrtho(left, right, bottom, top, near, far){
     // near*far > 0 !
     var rl1 = right-left;
     var rl2 = right+left;
     var tb1 = top-bottom;
     var tb2 = top+bottom;
     var fn1 = far-near;
     var fn2 = far+near;
     
     
     var result =  [ 2.0/rl1, 0,        0,         -rl2*1.0/rl1, 
                     0,       2.0/tb1,  0,         -1.0*tb2/tb1,  
                     0,       0,       -2.0/fn1,   -1.0*fn2/fn1,
                     0,       0,        0,          1];
     this.setByRows(result);
     
   }
   
   setFrustum(left, right, bottom, top, near, far){
     var rl1 = right-left;
     var rl2 = right+left;
     var tb1 = top-bottom;
     var tb2 = top+bottom;
     var fn1 = far-near;
     var fn2 = far+near;
     var  A = 1.0 * rl2 / rl1;    
     var  B = 1.0 * tb2 / tb1;    
     var  C = 1.0 * fn2 / fn1;    
     var D = 2.0 * far * near / fn1 ;
     var result =  [ 2.0*near/rl1, 0,            A, 0, 
                     0,            2.0*near/tb1, B, 0,  
                     0,            0,            C, D,
                     0,            0,            0, 1];
     this.setByRows(result);

   }

}


////////////////////////////////////////////////////
/// class VecN
////////////////////////////////////////////////////

class VecN {
  constructor(dim){
    this.dim = dim;
    this.array = new Array(dim);
    for(var i=0;i<dim;i++)
       this.array[i]=0.0;
  }
  
  getCopy(){
    var v = new VecN(this.dim);
    v.setByArray(this.array);
    return v;
  }
  
  setByArray(arr){
    
    this.array = arr.slice();

  }
  
  norm(){
    
    return Math.sqrt(VecN.dot(this,this));
    
  }
  
  normalize(){
    
    var s = this.norm();
    if (s!==0.0)
      for (var i=0;i<this.dim;i++)
         this.array[i]/=s;
    
  }
  
  ///////////////////////////////////////////////////////
  
  static add(vA,vB){
    var vC = new VecN(vA.dim);
    for(var i=0;i<vA.dim;i++)
       vC.array[i]=vA.array[i]+vB.array[i];
    return  vC;
  }
  
  static mMul(vM,vV){
    var vC = new VecN(vV.dim);
    for(var i=0;i<vV.dim;i++){
       var s = 0;
       for(var j=0;j<vV.dim;j++){
          s += vM.getAt(i,j)*vV.array[j];
       }
       vC.array[i] = s;
    }   
    return  vC;
  }
  
  static dot(vV1,vV2){
    
    var s = 0.0;
    for(var i=0;i<vV1.dim;i++)
       s += vV1.array[i]*vV2.array[i];
    return s;
  }

  
}

////////////////////////////////////////////////////
/// class Vec4
////////////////////////////////////////////////////

class Vec4 extends VecN {
  constructor(a,b,c,d){
    super(4);
    this.array[0] = a;
    this.array[1] = b;
    this.array[2] = c;
    this.array[3] = d;
  }
  setByPosition(x,y,z){
    this.array[0] = x;
    this.array[1] = y;
    this.array[2] = z;
    this.array[3] = 1;
  }
  setByDirection(x,y,z){
    this.array[0] = x;
    this.array[1] = y;
    this.array[2] = z;
    this.array[3] = 0;
  }
}

////////////////////////////////////////////////////
/// class Vec3
////////////////////////////////////////////////////


class Vec3 extends VecN {
  constructor(a,b,c){
    super(3);
    this.array[0] = a;
    this.array[1] = b;
    this.array[2] = c;
  }
  static cross(v1,v2){
    var a1 = v1.array;
    var a2 = v2.array;
    return new Vec3(a1[1] * a2[2] - a2[1]*a1[2],
                    -a1[0] * a2[2] + a2[0]*a1[2],
                    a1[0] * a2[1] - a2[0]*a1[1]);
  }
}

////////////////////////////////////////////////////
/// class Quaternion
////////////////////////////////////////////////////


class Quaternion{
  constructor(a=0,b=0,c=0,d=0){
    this.array = [a,b,c,d];
  }
  
  get a(){
    return this.array[0];
  }
  
  get b(){
    return this.array[1];
  }
  
  get c(){
    return this.array[2];
  }
  
  get d(){
    return this.array[3];
  }
  
  toString(){
    //console.log(this.array);
    var vs = ["*i","*j","*k"];
    var res = ""
    if (this.a!==0) res = res + this.a;
    var sn = "+";
    for (var i=1;i<4;i++){
      if (this.array[i]<0.0)  sn = "";
      else sn = "+";
      if (this.array[i]!==0.0) res = res + sn + this.array[i]+vs[i-1];
    }
    return  res;
    
  }
  
  setByArray(arr){
    this.array = arr.slice();
  }
  
  static add(qA,qB){
    var sum = new Quaternion(0,0,0,0);
    for(var i=0;i<4;i++)
       sum.array[i] = qA.array[i]+qB.array[i];
    return sum;
  }
  
  static qMul(qA,qB){
    var a1 = qA.array[0];var a2 = qB.array[0];
    var b1 = qA.array[1];var b2 = qB.array[1];
    var c1 = qA.array[2];var c2 = qB.array[2];
    var d1 = qA.array[3];var d2 = qB.array[3];
    
    var a = a1 * a2 - b1 * b2 - c1 * c2 - d1 * d2;
    var b = a1 * b2 + b1 * a2 + c1 * d2 - d1 * c2;
    var c = a1 * c2 - b1 * d2 + c1 * a2 + d1 * b2;
    var d = a1 * d2 + b1 * c2 - c1 * b2 + d1 * a2;
    
    
    return new Quaternion(a,b,c,d);
    
    
  }
  
  static qDiv(qA,qB){
    var cB = Quaternion.getConj(qB);
    var ls = qB.a * qB.a + qB.b*qB.b+qB.c*qB.c+qB.d*qB.d;
    return Quaternion.sMul((1.0/ls),cB);
    
  } 
  
  static sMul(s,qA){
    return new Quaternion(s*qA.a, s*qA.b, s*qA.c, s*qA.d);
  }

  static getConj(qA){
    return new Quaternion(qA.a, -qA.b, -qA.c, -qA.d);
  }
  
  static norm(qA){
    return Math.sqrt(qA.a*qA.a+qA.b*qA.b+qA.c*qA.c+qA.d*qA.d);
  }
  
  normalize(){
    var l = Quaternion.norm(this)*1.0;
    if (l!==0){
       this.a/=l;  this.b/=l; this.c/=l; this.d/=l;
    }
   
  }

  
}
/*
var m3 = new Mat3x3(0);
m3.setByRows([1,2,3,4,5,0,6,0,0]);

var m5 = new MatNxN(5,0);
m5.setByRows([10, 2, 3, 4, 8,11, 20, 3, 5, 9,12, 2, 40, 6, 10,13, 6, 5, 70, 11,14, 6, 6, 8, 12]);
 m5.setByRows([1, 0, 0, 0, 0, 
               0, 1, 0, 0, 0,
               0, 0, 2, 0, 0,
               0, 0, 0, 4, 0,
               0, 0, 0, 0, 8]);

*/


/*
var q2 = new Mat4x4(0);
q2.setOrtho(-10,10,-10,10,1,10);
console.log(VecN.mMul(q2, new Vec4(0,0,-10,1)));



console.log(q2.toString());

*/
/*
var vv = new Vec3(1,1,0);
var vv2 = new Vec3(0,1,1);
var vv3 = new Vec3(0,0,1);

var model = new Mat4x4(0);
model.setRotation(45,0,1,1);
console.log(model.toString());

console.log(MatNxN.getInverse(model).toString());
*/
/*
console.log(MatNxN.getDet(m5));
var invM5 = MatNxN.getInverse(m5);
console.log(invM5.toString());
console.log(MatNxN.mMul(invM5,m5).toString());

console.log(m3.det());

console.log(Vec3.cross(vv,vv2));
*/
