#include <iostream>
#include <cassert>
typedef long long ll;
using namespace std;

const static int P = 998244353;
const static int MAX = 900001;

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

// TODO change (from, to) to (p, p)

/**
 * @brief forward fast Fourier transform
 *
 * @note Both 'from' and 'to' are assumed to be have exactly MAX size.
 * @param n length of valid index for both from and to
 */
void *fft(const Mod *from, Mod *to, int n);

/**
 * @brief inverse fast Fourier transform
 *
 * @note Both 'from' and 'to' are assumed to be have exactly MAX size.
 * @param n length of valid index for both from and to
 */
void *ifft(const Mod *from, Mod *to, int n);

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
  Mod p[MAX]; // use when x is odd
  Mod q[MAX]; // use when x is else

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

  // TODO fft test

  // TODO ifft test

  // times2 test
  {
    Mod from[MAX], to[MAX];
    for(int i=0; i<5; i++){ from[i] = Mod(i+1); }
    for(int i=5; i<MAX; i++){ from[i] = Mod(0); }

    times2(from, to, 5);
    for(int j=0; j<20; j++)
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

