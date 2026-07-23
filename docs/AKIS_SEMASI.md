# CoreForge — Akış Şemaları

Bu doküman, kaynak kodun anlaşılmasını kolaylaştırmak amacıyla sistemin
ve ana hesap döngülerinin akış şemalarını içerir (Teknofest NDK
dokümantasyon şartı). Şemalar Mermaid formatındadır; GitHub üzerinde
doğrudan görüntülenir.

## 1 · Sistem mimarisi

```mermaid
flowchart TB
    subgraph UI["Streamlit Web Arayüzü (app.py)"]
        B[📚 Benchmark kütüphanesi]
        FD[🧬 Yakıt tasarımcısı<br/>zenginlik→kesit]
        M[🧪 Malzemeler<br/>makroskopik kesitler]
        CB[🧱 2D/3D Kor kurucu<br/>zonlar, yükleme deseni]
        S[⚡ Çöz & Sonuçlar]
        TH[🌡 Isıl-hidrolik]
        BU[🔥 Yakıt çevrimi<br/>+ çok-çevrim]
        TR[⏱ Transient + Xenon]
        T[🔧 Fizik araçları]
    end
    subgraph PY["Python katmanı"]
        R[runner.py<br/>girdi yazımı + çıktı ayrıştırma]
        X[xslib.py<br/>pin-hücre kesit üreteci +<br/>ters çözücü]
        BUP[burnup.py<br/>Bateman zinciri + letdown<br/>+ batch reload]
        K[kinetics.py<br/>nokta kinetiği + inhour<br/>+ Xe-135 transient]
        THM[thermal.py<br/>kapalı-kanal T-H]
        REP[report.py<br/>HTML rapor]
    end
    subgraph F["Fortran motoru (solver/coreforge.f90)"]
        E["2D/3D çok-gruplu difüzyon<br/>güç iterasyonu + red-black SOR<br/>OpenMP paralel"]
    end
    FD --> M
    M --> CB --> S
    S --> R --> E --> R --> S
    BU --> BUP --> R
    BUP --> X
    TR --> K
    TH --> THM
    S --> THM
    S --> REP
    B --> M
```

## 2 · Fortran özdeğer çözücüsü (bir k_eff çözümü)

```mermaid
flowchart TB
    A[girdi dosyasını oku<br/>NG,NX,NY,NZ, BC, malzemeler, MAP] --> V[doğrulama<br/>id aralığı, D>0, χ toplamı]
    V --> C["7-nokta katsayıları ÖN-HESAPLA<br/>harmonik ortalama D, Robin vakum"]
    C --> I["başlangıç: φ=1, F=νΣf·φ, normalize"]
    I --> O{dış iterasyon}
    O --> G[grup döngüsü g=1..NG]
    G --> SRC["kaynak: χ_g·F/k + saçılma"]
    SRC --> SOR["NINNER × red-black SOR süpürmesi<br/>(OpenMP, renk = i+j+k paritesi)"]
    SOR --> G
    G --> F2["F yenile, k = k·ΣF_yeni"]
    F2 --> N[normalize, artıklar dk, dS]
    N --> Q{dk<1e-7 ve dS<1e-5?}
    Q -- hayır --> O
    Q -- evet --> OUT["KEFF, F_q, nötron dengesi,<br/>flux.csv / power.csv"]
```

## 3 · Yakıt çevrimi (quasi-statik burnup)

```mermaid
flowchart TB
    A[designer yakıtlar → blok başına<br/>nüklit vektörü N] --> S1[akı çöz - Fortran]
    S1 --> XE["Xe-135/Sm-149 dengesi<br/>(mutlak akı: özgül güçten)"]
    XE --> LD{letdown?}
    LD -- evet --> CB["kritik boron ara<br/>(bisection, k=1)"]
    CB --> S2[akı çöz]
    LD -- hayır --> S2
    S2 --> EOC{çevrim sonu?<br/>otomatik: boron=0 ve k<1}
    EOC -- evet --> R["rapor: çevrim yanması, EFPD,<br/>yanma haritası, letdown eğrisi"]
    EOC -- hayır --> DEP["Bateman zinciri RK4:<br/>U-235↓, U-238→Pu-239→Pu-240→Pu-241,<br/>FP birikimi (blok blok)"]
    DEP --> XS["kesitleri N'den yeniden kur<br/>(xslib.cell_from_N)"]
    XS --> S1
```

## 4 · Transient (nokta kinetiği)

```mermaid
flowchart TB
    A["senaryo: step / rampa / scram<br/>ρ_dış(t), P₀, Λ"] --> B["y = [n, c₁..c₆] denge başlangıcı"]
    B --> L{zaman adımı}
    L --> R["ρ(t) = ρ_dış(t) + α_D·(T−T₀)"]
    R --> E["y ← expm(A(ρ)·dt)·y<br/>(donmuş-ρ KESİN çözüm)"]
    E --> T["yakıt sıcaklığı:<br/>M·cp dT/dt = (P−P₀) − M·cp(T−T₀)/τ_c"]
    T --> L
    L -- bitti --> OUT["P(t), ρ(t), T(t), tepe güç,<br/>asimptotik periyot ⇄ inhour"]
```

## 5 · Doğrulama zinciri (verify.py — 27 kontrol)

```mermaid
flowchart LR
    AN["Analitik<br/>k∞, çıplak kare,<br/>çıplak küp, inhour"] --> V((verify.py))
    BM["Uluslararası benchmark<br/>IAEA-2D, IAEA-3D,<br/>C5G7"] --> V
    INV["Değişmezlik<br/>2D ≡ extrude-3D"] --> V
    TR2["Fizik trendleri<br/>k∞(e), boron değeri,<br/>yanma limiti, Xe değeri"] --> V
    V --> P[PASS / FAIL raporu]
```
