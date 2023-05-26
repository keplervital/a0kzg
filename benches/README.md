# Domain Record Verification

The following explores the different alternatives for performing records verification with the newly proposed Name System.

## Benchmarks

The benchmarks shown are using `Apple M1 Pro CPU` and aim to evaluate the feasibility of a potential verkle tree implementation to perform verification of the domain records.

## Commit & Proof generation/verification benchmarks with KZG Affine Untrusted

The following table refers to times for the calculation of one commitment based on the set of values given. 

| Values |  Commit   |   Prove 1 |   Verify 1 |
|--------|-----------|-----------|------------|
|   32   | 3.9553 ms | 3.4846 ms |  3.8206 ms |
|   64   | 8.2745 ms | 7.1550 ms |  3.8176 ms | 
|  128   | 18.022 ms | 14.728 ms |  3.8245 ms | 
|  256   | 39.898 ms | 29.653 ms |  3.6923 ms | 
|  512   | 100.77 ms | 63.639 ms |  3.6994 ms | 

In the context of verkle trees the number of values could be interpreted as the width of the tree, the depth of the tree grows logarithmically based on the width (e.g 10000 elements in a tree with a widtyh of 32 would have a depth of `log32(10000)` which is 3). Larger widths of the tree result in a smaller proof size due to it's depth being smaller.

The following table estimates the size of the tree based on the number of elements and it's width:

| Elements |  Width | Height |     Insert 1 |     Insert All | Insert Elem/s | Proof Size |
|----------|--------|--------|--------------|----------------|---------------|------------|
| 1000000  | 32     |      4 |  115.4625 ms |   32.073 hours |     8.66      |  248 bytes |
| 1000000  | 64     |      4 |  466.1945 ms |       5.4 days |     2.15      |  248 bytes |
| 1000000  | 128    |      3 |  1903.206 ms |     22.03 days |     0.53      |  200 bytes |
| 1000000  | 256    |      3 |  7631.066 ms |     88.32 days |     0.13      |  200 bytes |
| 1000000  | 512    |      3 | 32683.938 ms |    378.29 days |     0.03      |  200 bytes |

**Insert 1 formula:** `(commit_ms + prove_1_ms * width)`

**KZG commitments:** Assuming `48-byte` (the commitments add up based on the depth of the tree)

**Merkle tree proof size:** Assuming `32-byte` hashes a merkle tree would have a proof size of around `750 bytes` with 100000 elements due to it's bigger depth

## Verkle tree vs Merkle tree 

Having a tree with the width of 32 seems to be a good balance between ahead of time computation and proof size, however, in the case of domain lookups using a merkle tree would allow for faster insert operations which would facilitate domain records management and the 750 bytes per proof for the example used would be generally only happening for the TLD operator canister which is less of an issue since that would be only one of the requests made by the custom domain resolver.

|        Scheme      |  Construct |    Update   | Proof Size | 
|--------------------|------------|-------------|------------|
| Merkle Tree        |    O(n)    |   O(log2 n) |  O(log2 n) |
| k-ary Verkle Tree  |    O(kn)   | O(k logk n) |  O(logk n) |

## DNSSEC

This approach relies on [RRSIG](https://datatracker.ietf.org/doc/html/rfc4034#section-3) records and [DNSKEY](https://datatracker.ietf.org/doc/html/rfc4034#section-2) records to establish a chain of trust across the domain hierarchy with a trusted anchor that is hard coded in the resolver. For the resolver to be able to check the signatures available in the RRSIG records of the record type it requested for, it needs the public key that is available in the DNSKEY record, those records are returned as extra information of the domain record query and the size varies based on how many public keys are associated with the domain.  

Latest DNSSEC schemes uses ECDSA keys with the P-256 curve, those keys use `32 bytes` and generated signatures thgat are `64 bytes` long that when base64 encoded increase the size to around `86 bytes`. Name lookups generally associate only one record type per request which means that when DNSSEC is enable N DNSKEY records and the RRSIG for the requested record type is returned which translates to around `118 bytes`.

Verification times with ECDSA are also expected to be constant which would be fast when compared with merkle or verkle trees that increase with the depth of the tree.

### On-chain ECDSA signatures 

Some Internet Computer subnets already allow for signatures to be created from the canister in a secure manner, however certain limitations apply:

- The is an average cost of 0.03$ per signature
- Throughtput is currently around 1 signature per second per subnet
- The curve used is not support by DNSSEC, we would need to support P-256
- ECDSA signatures are not available to all subnets

Those limitations are most troublesome in the case of DAOs that would require signatures to be performed on chain.

## Conclusion

DNSSEC is the preferred solution, not only it would be compatible with Web 2 it would also result in less data transfered per name lookup and shorter verification times on larger datasets.

To enable it two main limitations need to be resolved:

- IC adds support to the P-256 curve
- Enable ECDSA signatures to all subnets

The envisioned solution would enable DAOs to perform signatures on chain and for centrally owned canisters the controllers can perform the signatures off-chain and add the records in the same manner as other records are updated.

When DNSSEC is performed through the custom domain resolver the Domain Root canister public key is the hardcoded trust anchor.

**Note:** If those limitations are not resolved then a Merkle tree seems to be the custom solution that is not supported by Web 2 but would provide a good balance between data transfer size of name lookups and tree generation times as opposed to a Verkle tree with `width=32`.

## References

[Verkle Trees by John Kuszmaul](https://math.mit.edu/research/highschool/primes/materials/2018/Kuszmaul.pdf)

[Verkle Trees: Ver(y Short Mer)kle Trees by John Kuszmaul](https://math.mit.edu/research/highschool/primes/materials/2019/conf/12-5-Kuszmaul.pdf)

[Verkle trees by Vitalik Buterin](https://vitalik.ca/general/2021/06/18/verkle.html)
