#include <iostream>
#include <cassert>
#include <utility>
typedef long long ll;
using namespace std;

const static int P = 998244353;
const static int MAX = 1 << 20;

// class for Modulo P operation
struct Mod
{
  Mod()
  {
    // Do nothing
  }
  Mod(int val)
  {
    _val = val % P;
  }

  Mod(ll val)
  {
    _val = (int) ( val % (ll)P );
  }

  int val(){ return _val; }

  Mod operator+(Mod rhs)
  {
    int temp = this->val() + rhs.val();
    temp %= P;
    return Mod(temp);
  }

  Mod operator-(Mod rhs)
  {
    int temp = this->val() - rhs.val() + P;
    temp %= P;
    return Mod(temp);
  }

  Mod operator*(Mod rhs)
  {
    ll temp = ( (ll)(this->val()) * (ll)(rhs.val()) ) % (ll)P;
    return Mod(temp);
  }

  Mod operator/(Mod rhs)
  {
    return (*this) * inverse(rhs);
  }

  bool operator==(Mod rhs)
  {
    return *this == rhs.val();
  }

  bool operator==(int rhs)
  {
    return (this->val() == rhs) ? true : false;
  }

  static Mod power(Mod arg, int pow)
  {
    if(pow == 1){ return arg; }
    if(pow % 2)
    {
      Mod temp = power(arg, pow/2);
      return temp * temp * arg;
    }
    else
    {
      Mod temp = power(arg, pow/2);
      return temp * temp;
    }
  }

  private:
  Mod inverse(Mod arg)
  {
    assert( arg.val() != 0 );
    return power(arg, P-2);
  }

  int _val;
};

const static Mod G(3);

/**
 * @brief generate rev according to length n
 *
 * So called Butterfly exchange. It is reordering of [0,n)
 * In case of n = 8, rev is 04261537
 * In case of n = 16, rev is 0 8 4 12 2 10 6 14 1 9 5 13 3 11 7 15
 *
 * @param n Length of rev to make. It SHOULD be power of 2.
 */
void make_rev(int *rev, int n)
{
  int len = 0;
  for(int m=n; !(m & 1); m >>= 1)
    len++;
  rev[0] = 0;
  for (int i = 1; i < n; i++)
    rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (len - 1));
}

// Mod version of number_theoretic_transform()
void NTT(Mod *p, int n, bool is_inverse = false)
{
  int rev[MAX];
  make_rev(rev, n);
  for (int i = 0; i < n; i++)
    if (i < rev[i])
      std::swap(p[i], p[rev[i]]);
  for (int j = 1; j < n; j <<= 1)
  {
    Mod wn1, w, t0, t1;
    wn1 = Mod::power(G, (P - 1) / (j << 1));
    if (is_inverse)
      wn1 = Mod::power(wn1, P - 2);
    for (int i = 0; i < n; i += j << 1)
    {
      w = Mod(1);
      for (int k = 0; k < j; k++)
      {
        t0 = p[i + k];
        t1 = w * p[i + j + k];
        p[i + k] = t0 + t1;
        p[i + j + k] = t0 - t1;
        w = w * wn1;
      }
    }
  }
  if (is_inverse)
  {
    for (int i = 0; i < n; i++)
      p[i] = p[i] / Mod(n);
  }
}

/**
 * @brief generate from P[,x] to P[,2x]
 *
 * It basically does convolution using NTT.
 * @note Both 'from' and 'to' are assumed to be have exactly MAX size.
 * @param flen length of valid (non-zero) 'from' index
 */
void times2(const Mod *from, Mod *to, int flen);

/**
 * @brief generate from P[,x] to P[,x+1]
 *
 * It basically implements recurrent relation.
 * @note Both 'from' and 'to' are assumed to be have exactly MAX size.
 * @param flen length of valid (non-zero) 'from' index
 * @param used list of used digits, e.g. 4 0 2
 */
void next(const Mod *from, Mod *to, int flen, int *used);

int prob()
{
  // input
  int n, m, k; cin >> n >> k;
  m = n / 2;
  int d[10];
  int dmax = -1;
  for (int i=0; i<k; i++)
  {
    cin >> d[i];
    if (dmax < d[i]) { dmax = d[i]; }
  }

  int x = 1;
  // p[], q[] are # of l-digit s.t. sum is x
  static Mod p[MAX]; // use when x is odd TODO delete static?
  static Mod q[MAX]; // use when x is else

  // initialize
  for(int s=0; s<=dmax; s++){ p[s] = Mod(0); }
  for(int i=0; i<k; i++){ p[d[i]] = Mod(1); }

  Mod *from;
  Mod *to;
  while(x != m)
  {
    x++;
    if(x % 2)
    {
      from = q;
      to = p;
    }
    else
    {
      from = p;
      to = q;
    }

    for(int s=0; s<= x*dmax; s++)
    {
      to[s] = Mod(0);
      for(int i=0; i<k; i++)
      {
        if(s - d[i] >= 0 && s - d[i] <= (x-1)*dmax)
        {
          to[s] = to[s] + from[s - d[i]];
        }
      }
    }
  }

  assert( x == m );
  Mod count(0);
  Mod *head;

  if (x % 2)
  {
    head = p;
  }
  else
  {
    head = q;
  }

  for (int s=0; s<= x*dmax; s++)
  {
    count = count + head[s] * head[s];
  }

  int ans = count.val();

  // solve prob

  cout << ans << endl;
  return ans;
}

void test()
{
  // Mod test
  {
    assert( Mod(P-1) + Mod(P-1) == P-2 );
    assert( Mod(1) - Mod(2) == P-1 );
    assert( Mod(0) - Mod(P-1) == 1 );
    assert( Mod(P-1) * Mod(P-1) == 1 );
    assert( Mod(1) + Mod(P-1) * Mod(P-1) == 2 );
    assert( Mod(P-1) * Mod(P-1) * Mod(P-1) == P-1 );
    assert( Mod(1) / Mod(P-1) == P-1 );

    assert( Mod(7) / Mod(8) == 124780545);
    assert( Mod(1)/Mod(2) + Mod(1)/Mod(3) == Mod(5)/Mod(6) );
    assert( Mod(1)/Mod(2) - Mod(1)/Mod(3) == Mod(1)/Mod(6) );
    assert( (Mod(1)/Mod(2)) / Mod(3) == Mod(1)/Mod(6) );
  }

  // make_rev test
  {
    int n = 16;
    int r[16];
    make_rev(r, n);

    int expected[16] = {0,8,4,12,2,10,6,14,1,9,5,13,3,11,7,15};
    for (int i=0; i<n; i++)
    {
      assert( r[i] == expected[i] );
    }
  }

  // NTT test
  {
    int n = 16;
    Mod p[16];
    for(int i=0; i<n; i++) { p[i] = Mod(i); }
    NTT(p, n);

    int expected[16] =
    {
      120, 16886715, 790357655, 115058691, 692669736, 306777988, 403262520, 432660095, 998244345, 565584242, 594981817, 691466349, 305574601, 883185646, 207886682, 981357622
    };
    for(int i=0; i<n; i++)
    {
      assert( p[i] == expected[i] );
    }
  }

  // times2 test
  {
    int n = 5;
    Mod from[MAX], to[MAX]; // TODO static?
    for(int i=0; i<n; i++){ from[i] = Mod(i); }

    times2(from, to, n);
    for(int j=0; j<2*n-1; j++)
    {
      Mod conv(0);
      for(int i=0; i<=j; i++)
      {
        conv = conv + from[i] * from[j-i];
      }
      assert( to[j] ==  conv );
    }
  }

  // TODO next test

  // self-judge by input.txt
  {
    freopen("input.txt", "r", stdin);

    int T; cin >> T;
    for (int tc=1; tc<=T; tc++)
    {
      cout << "test " << tc << endl;
      int ans = prob();
      int expected; cin >> expected;

      assert( ans == expected );
    }
  }
}

int main()
{
  test();
  //prob();

  return 0;
}

